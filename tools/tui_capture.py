#!/usr/bin/env python3
"""Drive rfam_top inside a real pty and capture one rendered frame.

Test harness only: it exists so the curses layout can be inspected without an
interactive terminal (CI, agent sessions, this repository's WSL image).
"""
import fcntl
import os
import pty
import select
import struct
import sys
import termios
import time

ROWS, COLS = 44, 140


def main() -> int:
    argv = sys.argv[1:]
    if not argv:
        print("usage: tui_capture.py <program> [args...]")
        return 2
    run_seconds = float(os.environ.get("CAPTURE_SECONDS", "4"))

    pid, fd = pty.fork()
    if pid == 0:
        os.execvp(argv[0], argv)

    fcntl.ioctl(fd, termios.TIOCSWINSZ, struct.pack("HHHH", ROWS, COLS, 0, 0))
    os.kill(pid, 28)  # SIGWINCH so ncurses picks up the size

    chunks = []
    deadline = time.time() + run_seconds
    while time.time() < deadline:
        ready, _, _ = select.select([fd], [], [], 0.2)
        if not ready:
            continue
        try:
            data = os.read(fd, 65536)
        except OSError:
            break
        if not data:
            break
        chunks.append(data)

    try:
        os.write(fd, b"q")
    except OSError:
        pass
    time.sleep(0.4)
    try:
        while True:
            ready, _, _ = select.select([fd], [], [], 0.2)
            if not ready:
                break
            data = os.read(fd, 65536)
            if not data:
                break
            chunks.append(data)
    except OSError:
        pass

    try:
        os.waitpid(pid, os.WNOHANG)
    except ChildProcessError:
        pass

    sys.stdout.buffer.write(b"".join(chunks))
    return 0


if __name__ == "__main__":
    sys.exit(main())
