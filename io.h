#ifndef IO_H
#define IO_H

/* types of user command-line input */
typedef enum {
  INT,
  DOUBLE,
  STR,
  NA
} ARG_TYPE;

#endif

int findarg(const char *argname, ARG_TYPE type, void *val, int argc, char **argv);

// This part is copied from http://www-users.cs.umn.edu/~saad/software/EVSL/
