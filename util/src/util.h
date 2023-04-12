/*
 *  Copyright (C) 2019 Tyson B. Littenberg (MSFC-ST12), Neil J. Cornish
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/*
 * This is a set of convenience functions to make certain patterns that gcc -Wall
 * does not like have well defined behavior. They generally avoid the problems gcc -Wall is
 * warning about by reporting it loudly at runtime.
 * 
 * Code that wants to handle these file and string operations gracefully should NOT use these
 * functions, and ought carefully handle the returns and conditions on the underlying library 
 * functions under -Wall.
 * 
 * Most of these have the happy-path implemented by a macro. This is intentional both to allow
 * The compiler to default to inlining, and easily optimize out any checks it can statically 
 * prove will never hit. Macros also allow the easy use of __FILE__ , __LINE__, and __func__ 
 * for usable error messages. These are important because the exit pathway in panic() is not 
 * able to print a stack trace, and by design these macros should be all over the code.
 * 
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

/* 
Preprocessor gadget to allow calls to variadic macros without any arguments. 
 - It expands to comma followed by whatever args you pass in.
 - When you pass no arguments the token concatenation can't so the whole thing expands to nothing

IN C++20 there is __VA_OPT__ to do this, but we want to work on compilers that don't have that yet
*/
#define VA_ARGS(...) , ##__VA_ARGS__


/**
 * Unsafe file reading functions
 * 
 * Replacements for fscanf/fread that panics if the file doesn't scan
 * It eats the return value though, so if you need that, don't use these. 
 * 
 * This converts the condition -Wunused-result is worried about into a helpful runtime error
*/
#define ufscanf(stream, format, ...) ({\
    int retval = fscanf(stream, format VA_ARGS(__VA_ARGS__));\
    if(retval < 0) panic("ufscanf error: fscanf returned %d from __FILE__:__LINE__ in __func__(). Exiting now!", retval);\
})

#define ufread(ptr, size, nmemb, stream) ({\
    size_t retval = fread(ptr, size, nmemb, stream);\
    if(retval != nmemb) panic("ufread error: fread didnt read what caller expected. Returned %d from __FILE__:__LINE__ in __func__(). Exiting now!", retval);\
})


/**
 * Path manipulation functions
*/

/* 
 * Replacement for MAX_PATH inheriting its many problems
 *
 * It's important to note that even PATH_MAX in limits.h is not really the maximum length for the path
 * because it is f/s dependent, only available at runtime, and depends on the working directory, because
 * maximum path lengths are (in most f/s) limitations on the length of an absolute path. 
 */
#define PATH_BUFSIZE 1024

/**
 * Works like sprintf, but insists caller pass a buffer of PATH_BUFSIZE. 
 * Panics if anything goes wrong. No return value.
 * 
 * This construction avoids -Wformat-overflow and -Wformat-truncation messages from the compiler 
 * by turning both conditions into runtime errors.
 *   
 */
#define pathprintf(buffer, format, ...) ({\
    int retval = snprintf(0, 0, format VA_ARGS(__VA_ARGS__));\
    if(retval < 0 || retval >= PATH_BUFSIZE) {\
        panic("pathprintf error: snprintf returned %d from __FILE__:__LINE__ in __func__(). Exiting now!", retval);\
    }\
    retval = snprintf(buffer, PATH_BUFSIZE, format VA_ARGS(__VA_ARGS__));\
    if(retval < 0 || retval >= PATH_BUFSIZE) {\
        panic("pathprintf error: snprintf returned %d from __FILE__:__LINE__ in __func__(). Exiting now!", retval);\
    }\
})

/**
 * This construction is used to avoid -Wrestrict in situations where you want to sprintf(buf, "%s/<format stuff>", buf, <va_args>)
 * You instead call pathappendprintf(buf, "/<format stuff>", <va_args>)
 * 
*/
#define pathappendprintf(buffer, format, ...) ({\
    size_t orig_len = strnlen(buffer, PATH_BUFSIZE);\
    if(orig_len >= PATH_BUFSIZE) {\
        panic("pathappendprintf error: Trying to append to a path that is already PATH_BUFSIZE or greater. Called from  __FILE__:__LINE__ in __func__(). Exiting now!");\
    }\
    int retval = snprintf(0, 0, format VA_ARGS(__VA_ARGS__));\
    if(retval < 0 || retval >= PATH_BUFSIZE - orig_len) {\
        panic("pathappendprintf error: snprintf to null returned %d from __FILE__:__LINE__ in __func__(). Exiting now!", retval);\
    }\
    char tmpbuf[PATH_BUFSIZE];\
    strncpy(tmpbuf, buffer, orig_len);\
    retval = snprintf(buffer, PATH_BUFSIZE, "%s"format, tmpbuf VA_ARGS(__VA_ARGS__));\
    if(retval < 0 || retval >= PATH_BUFSIZE) {\
        panic("pathappendprintf error: snprintf returned %d from __FILE__:__LINE__ in __func__(). Exiting now!", retval);\
    }\
})

/**
 * Used by many of these macros to generate a runtime error to stderr
 * This function does not return.
 */
void panic(char* message, ...);


