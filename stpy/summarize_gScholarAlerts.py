"""
This scripts is aiming at automating the process of summarizing the publications from Google Scholar Alerts.

It connects to a Gmail account, fetches all/unread emails from the "gScholarAlerts" label for the last 7 days, extracts publication details,
and saves the summary to a CSV file. It also retrieves additional information from the DOI and URL of the publications.
The script is designed to work with Google Scholar Alerts, which send notifications about new publications based on user-defined search queries.

The script uses the following libraries:
- imaplib: For connecting to the Gmail IMAP server and fetching emails.
- email: For parsing the email content.
- pandas: For data manipulation and saving to CSV.
- BeautifulSoup: For parsing HTML content.
- requests: For making HTTP requests to fetch publication details from DOI and URL.
- fitz (PyMuPDF): For extracting text from PDF files.
The script is designed to be run as a standalone program, and it requires the following environment variables to be set:
- GMAIL_USERNAME: The Gmail username (email address).
- GMAIL_PASSWORD: The Gmail password (or app password if 2FA is enabled).
It is recommended to use an app password for security reasons if 2FA is enabled on the Gmail account.

"""

import email
import imaplib
import logging
import os
import re
from datetime import datetime, timedelta
from io import BytesIO

import fitz  # PyMuPDF
import pandas as pd
import requests
from bs4 import BeautifulSoup

# import openpyxl


# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def fetch_unread_emails(username, password, since_days, label="gScholarAlerts"):
    # label:gscholaralerts
    try:
        # Connect to Gmail IMAP server
        mail = imaplib.IMAP4_SSL("imap.gmail.com")
        mail.login(username, password)

        # Select emails from the specific label (default: 'gscholaralerts')
        status, _ = mail.select(f'"{label}"')  # Ensure label exists & use double quotes

        if status != "OK":
            logging.error(f"Failed to select label: {label}")
            return None, []

        # Calculate date filter (IMAP format: DD-MMM-YYYY)
        since_date = (datetime.today() - timedelta(days=since_days)).strftime(
            "%d-%b-%Y"
        )

        # Search query within the label instead of inbox
        # search_query = f'(UNSEEN FROM "scholaralerts-noreply@google.com" SINCE {since_date})'
        search_query = f'(ALL FROM "scholaralerts-noreply@google.com" SINCE {since_date})'  # change "ALL" to "UNSEEN" to get all unseen emails
        status, messages = mail.search(None, search_query)

        if status != "OK":
            logging.error(f"Search command failed for label: {label}")
            return None, []

        email_ids = messages[0].split()
        logging.info(
            f"Fetched {len(email_ids)} emails from label '{label}' since {since_date}."
        )

        return mail, email_ids
    except Exception as e:
        logging.error(f"Failed to fetch emails: {e}")
        return None, []


def parse_email(mail, email_id):
    try:
        status, msg_data = mail.fetch(email_id, "(RFC822)")
        msg = email.message_from_bytes(msg_data[0][1])
        if msg.is_multipart():
            for part in msg.walk():
                if part.get_content_type() == "text/html":
                    return part.get_payload(decode=True)
        else:
            return msg.get_payload(decode=True)
    except Exception as e:
        logging.error(f"Failed to parse email {email_id}: {e}")
        return None


def extract_publication_details(html_content):
    try:
        soup = BeautifulSoup(html_content, "html.parser")
        publications = [
            {
                "href": a["href"],
                "title": a.get_text(),
                "author - source": a.find_parent("h3")
                .find_next_sibling("div")
                .get_text(strip=True),
            }
            for a in soup.find_all("a", class_="gse_alrt_title")
        ]
        logging.info(f"Extracted {len(publications)} publications.")
        return publications
    except Exception as e:
        logging.error(f"Failed to extract publication details: {e}")
        return pd.DataFrame()


def get_more_info(publications):
    # Retrieve other information from DOI
    publications_df = pd.DataFrame(publications)
    # drop duplicates
    # print(publications_df.columns)
    publications_df.drop_duplicates(subset=["href"], inplace=True)
    publications_df["url"] = publications_df["href"].apply(
        lambda x: x.split("scholar_url?url=")[1].split("&")[0]
    )

    doi_pattern = re.compile(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+", re.IGNORECASE)
    publications_df["doi"] = publications_df["href"].apply(
        lambda x: doi_pattern.search(x).group(0) if doi_pattern.search(x) else None
    )

    # Get more information from DOI
    publications_df["info_from_doi"] = publications_df["doi"].apply(
        lambda x: get_publication_details(x) if x else None
    )
    # Get more information from url
    df_expanded = publications_df.apply(
        lambda row: extract_publication_info(row["url"]), axis=1
    ).apply(pd.Series)
    # Combine with original DataFrame
    publications_df["info_from_url"] = publications_df["url"].apply(
        lambda x: extract_publication_info(x) if x else None
    )
    return publications_df


def get_publication_details(doi):
    """
    Retrieve publication details including abstract using DOI from the CrossRef API.
    """
    base_url = "https://api.crossref.org/works/"
    full_url = f"{base_url}{doi}"

    try:
        response = requests.get(full_url)
        response.raise_for_status()  # Raise an error if the request fails

        data = response.json()["message"]

        # Extract details
        title = data.get("title", ["No title available"])[0]
        authors = [
            author["given"] + " " + author["family"]
            for author in data.get("author", [])
        ]
        journal = data.get("container-title", ["No journal info"])[0]
        publication_year = data.get("published-print", {}).get("date-parts", [[None]])[
            0
        ][0]
        abstract = data.get("abstract", "No abstract available")  # Abstract field

        return {
            "DOI": doi,
            "Title": title,
            "Authors": authors,
            "Journal": journal,
            "Publication Year": publication_year,
            "Abstract": abstract,
        }

    except requests.exceptions.RequestException as e:
        print(f"Error retrieving DOI information: {e}")
        return None


# # Example usage
# doi = "10.1038/s41586-020-2649-2"  # Replace with any valid DOI
# publication_info = get_publication_details(doi)

# if publication_info:
#     print(publication_info)


def is_pdf_url(url):
    """
    Check if the URL points to a PDF file by examining headers.
    """
    response = requests.head(url, allow_redirects=True)
    return response.headers.get("Content-Type") == "application/pdf"


def extract_text_from_pdf(pdf_url):
    """
    Download a PDF from a URL and extract its text without saving it locally.
    """
    response = requests.get(pdf_url)
    if response.status_code == 200:
        pdf_stream = BytesIO(response.content)  # Load PDF into memory
        doc = fitz.open(stream=pdf_stream, filetype="pdf")  # Open PDF in memory
        text = "\n".join(
            [page.get_text("text") for page in doc]
        )  # Extract text from all pages
        doc.close()  # Close the document
        pdf_stream.close()
        return text
    else:
        print("Failed to fetch PDF.")
        return None


def extract_text_from_webpage(web_url):
    """
    Extract text from an HTML webpage.
    """
    headers = {"User-Agent": "Mozilla/5.0"}
    response = requests.get(web_url, headers=headers)

    if response.status_code != 200:
        print(f"Failed to fetch webpage: {web_url}")
        return None

    soup = BeautifulSoup(response.text, "html.parser")

    # Find the "Abstract" element
    abstract_tag = soup.find(
        lambda tag: tag.name in ["h2", "h3", "strong"]
        and "abstract" in tag.text.lower()
    )

    if abstract_tag:
        # Get the next sibling that contains text
        next_tag = abstract_tag.find_next_sibling()
        abstract_text = next_tag.text.strip() if next_tag else "No abstract found"
    else:
        abstract_text = "Abstract heading not found"

    return {"Abstract": abstract_text}


def extract_publication_info(url):
    """
    Extract publication details based on whether the URL is a PDF or a web page.
    """
    if is_pdf_url(url):
        text = extract_text_from_pdf(url)
        if text:
            # Extract title (assuming it's capitalized and at the start)
            title_match = re.search(r"(?:(?:[A-Z][a-z]+)\s+){3,}", text)
            title = title_match.group(0) if title_match else "Title not found"

            # Extract authors (assuming they appear in a format like: "Author 1, Author 2, ...")
            authors_match = re.search(r"(?:(?:[A-Za-z\-]+),\s*){2,}", text)
            authors = authors_match.group(0) if authors_match else "Authors not found"

            # Extract abstract (assuming it starts with 'Abstract' and ends before 'Introduction')
            # abstract_match = re.search(r"Abstract([\s\S]*?)Introduction", text, re.IGNORECASE)
            abstract_match = re.search(
                r"Abstract([\s\S]+?)(?:\n[A-Z][a-z]+\s*\n|\Z)", text, re.IGNORECASE
            )

            abstract = (
                abstract_match.group(1).strip()
                if abstract_match
                else "Abstract not found"
            )

            return {"Title": title, "Authors": authors, "Abstract": abstract}
        else:
            return None
    else:
        return extract_text_from_webpage(url)


# Example usage
# pdf_url = "https://example.com/sample.pdf"  # Replace with a valid PDF URL
# pdf_info = extract_text_from_webpage(pdf_url)


def summarize_publications(username, password, since_days):
    mail, email_ids = fetch_unread_emails(username, password, since_days)
    since_date = (datetime.today() - timedelta(days=since_days)).strftime("%d-%b-%Y")
    if not mail:
        logging.error("Failed to connect to mail server.")
        return
    all_publications = []
    for email_id in email_ids:
        html_content = parse_email(mail, email_id)
        if html_content:
            publications = extract_publication_details(html_content)
            all_publications += publications
    mail.logout()
    if all_publications:
        # Clean up the DataFrame
        # print(all_publications)
        df = get_more_info(all_publications)
        # Save the DataFrame to CSV and Excel
        df.to_csv(
            f'publications_summary{since_date}to{datetime.today().strftime("%d-%b-%Y")}.csv',
            index=False,
        )
        logging.info("Summary saved to publications_summary.csv")
    else:
        logging.info("No publications found.")


if __name__ == "__main__":
    username = os.getenv("GMAIL_USERNAME")
    password = os.getenv("GMAIL_PASSWORD")
    since_days = os.getenv("SINCE_DAYS", 7)  # Default to 7 days if not set
    if username and password:
        logging.info("Starting to summarize publications...")
        summarize_publications(username, password, since_days)
        logging.info("Summary completed.")
        # logging.info("No unread emails found.")
    else:
        logging.error(
            "GMAIL_USERNAME and GMAIL_PASSWORD environment variables not set."
        )
