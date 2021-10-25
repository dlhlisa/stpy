# Tmux

## [A Quick and Easy Guide to tmux](https://www.hamvocke.com/blog/a-quick-and-easy-guide-to-tmux/)

### Installation
- `sudo apt-get install tmux` (Ubuntu and derivatives) 
- `brew install tmux` (Mac)

### Basic usage
- `tmux` # start new session
- **`Ctrl + b` + %** # split panes
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

## [Source code](https://github.com/tmux/tmux)
