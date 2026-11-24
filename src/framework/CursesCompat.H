#ifndef CursesCompat_H
#define CursesCompat_H

// The dashboard needs ncurses. Two build modes exist because they are two
// different machines, not two different programs:
//
//   RFAM_HAVE_CURSES_HEADERS   system ncurses development headers are present
//                              (A406, and any normal build host)
//   RFAM_CURSES_SHIM           only the runtime soname exists (the WSL image
//                              this repository is developed on has
//                              libncursesw.so.6 but no headers and no sudo)
//
// In shim mode every symbol below is still resolved against the real ncurses
// runtime; only the declarations are supplied locally. WINDOW stays opaque and
// all access goes through functions, never through the WINDOW struct layout, so
// this does not pin an ncurses version.

#ifdef RFAM_HAVE_CURSES_HEADERS

#include <curses.h>

#else

#include <cstdarg>

extern "C"
{

typedef struct _win_st WINDOW;
typedef unsigned int chtype;
typedef unsigned int attr_t;

typedef short mmask_t;

extern WINDOW* stdscr;
extern int LINES;
extern int COLS;
extern int COLORS;
extern int COLOR_PAIRS;

WINDOW* initscr(void);
int endwin(void);
int isendwin(void);
int cbreak(void);
int noecho(void);
int nonl(void);
int keypad(WINDOW*, bool);
int nodelay(WINDOW*, bool);
int wtimeout(WINDOW*, int);
int curs_set(int);
int start_color(void);
int has_colors(void);
int assume_default_colors(int, int);
int init_pair(short, short, short);
int flushinp(void);

int wmove(WINDOW*, int, int);
int waddch(WINDOW*, const chtype);
int waddnstr(WINDOW*, const char*, int);
int wattron(WINDOW*, int);
int wattroff(WINDOW*, int);
int werase(WINDOW*);
int wclrtoeol(WINDOW*);
int wrefresh(WINDOW*);
int wnoutrefresh(WINDOW*);
int doupdate(void);
int wgetch(WINDOW*);
int getmaxy(WINDOW*);
int getmaxx(WINDOW*);
int resize_term(int, int);

}

// ncurses.h constants, from the stable ABI (NCURSES_ATTR_SHIFT == 8).
#define RFAM_A_NORMAL    0
#define RFAM_A_BOLD      (1U << 21)
#define RFAM_A_DIM       (1U << 20)
#define RFAM_A_REVERSE   (1U << 18)
#define RFAM_A_UNDERLINE (1U << 17)

#ifndef COLOR_BLACK
#define COLOR_BLACK   0
#define COLOR_RED     1
#define COLOR_GREEN   2
#define COLOR_YELLOW  3
#define COLOR_BLUE    4
#define COLOR_MAGENTA 5
#define COLOR_CYAN    6
#define COLOR_WHITE   7
#endif

#ifndef COLOR_PAIR
#define COLOR_PAIR(n) (static_cast<int>((n) << 8))
#endif
#ifndef A_BOLD
#define A_BOLD      RFAM_A_BOLD
#define A_DIM       RFAM_A_DIM
#define A_REVERSE   RFAM_A_REVERSE
#define A_NORMAL    RFAM_A_NORMAL
#define A_UNDERLINE RFAM_A_UNDERLINE
#endif

#ifndef KEY_UP
#define KEY_UP    0403
#define KEY_DOWN  0402
#define KEY_LEFT  0404
#define KEY_RIGHT 0405
#define KEY_NPAGE 0522
#define KEY_PPAGE 0523
#define KEY_RESIZE 0632
#define KEY_ENTER 0343
#define KEY_BACKSPACE 0407
#endif

#ifndef ERR
#define ERR (-1)
#endif

#endif // RFAM_HAVE_CURSES_HEADERS

#endif
