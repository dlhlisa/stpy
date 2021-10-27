# Tmux

## [A Quick and Easy Guide to tmux](https://www.hamvocke.com/blog/a-quick-and-easy-guide-to-tmux/)

### Installation
- `sudo apt-get install tmux` (Ubuntu and derivatives) 
- `brew install tmux` (Mac)

### Basic usage
- `tmux` # start new session
- **`Ctrl + b` + %** or **`Ctrl + b` + "** # split panes
- C-b <arrow key> switching to a different pane
- `exit` or `Ctrl + d` closing panes
- `C-b c` creating new windows  
- `C-b p` to switch to the previous window, or `C-b n` to switch to the next window
- `C-b d` to detach current session, or `C-b D` give you a choice which of your sessions you want to detach
- `tmux ls`
- `tmux attach -t 0` # -t 0 tells wich session to attach to
- `tmux new -s database` # create a new session with the name 'databast'
- `tmux rename-ression -t 0 database` # rename your existing session
- `tmux attach -t database` # attach named session
- `C-b ?` to see a list of all available commands and start experimenting
  - C-b z: make a pane go full screen. Hit C-b z again to shrink it back to its previous size
  - C-b C-<arrow key>: Resize pane in direction of <arrow key>4
  - C-b ,: Rename the current window

### [Source code](https://github.com/tmux/tmux)
  
# [Screen](https://www.gnu.org/software/screen/manual/screen.html

## [How To Use Linux Screen](https://linuxize.com/post/how-to-use-linux-screen/)
  
### Installation
- `sudo apt update`, `sudo apt install screen` (Ubuntu and derivatives) 
- `sudo yum install screen` (Mac)

### Basic usage
- Ctrl+a c Create a new window (with shell).
- Ctrl+a " List all windows.
- Ctrl+a 0 Switch to window 0 (by number).
- Ctrl+a A Rename the current window.
- Ctrl+a S Split current region horizontally into two regions.
- Ctrl+a | Split current region vertically into two regions.
- Ctrl+a tab Switch the input focus to the next region.
- Ctrl+a Ctrl+a Toggle between the current and previous windows
- Ctrl+a Q Close all regions but the current one.
- Ctrl+a X Close the current region.
- Ctrl+a d  detach from the screen session.
- `screen -ls` list all the screen sessions
- `screen -r ID/named_session` to resume your screen session.
- `screen` # start new session or `screen -S named_session` to start new named session.
- **`Ctrl + a` + |**  or **`Ctrl + b` + %** # split panes
- Ctrl + a + <tab key> switching to a different pane

## Comparison between Tmux and Screen

Tmux and Screen approximately serves a similar purpose. Both programs build a virtual Terminal inside a single Terminal, allows you to switch among virtual Terminals and let you attach and reattach the virtual Terminals when your network connection is disrupted. Both programs operate by building separate processes which they name differently.

However, there are some differences as well among these two tools. Tmux has a BSD license while the Screen has GNU GPL. Tmux is more user-friendly than the Screen and contains a nice status bar with some info in it. Tmux features automatic window renaming while the Screen lacks this feature. The Screen allows session sharing with other users while Tmux does not. That is the great feature that Tmux lacks.
  
[Reference can be found here](https://linuxhint.com/tmux_vs_screen/)


