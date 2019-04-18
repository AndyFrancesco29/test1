/*
written by Guido Perrone (perrone@polito.it)
*/

#define __gERR_H

typedef enum
    {
    zero,                         // no errors
    memory,                       // memory allocation
    invalid,                      // invalid operation
    few_points,                   // few points
    too_points,                   // too many points
    fail,                         // requested operation failed
    exists,                       // already exists
    eof,                          // end of file
    doubt,                        // is it an error ?
    out_range,                    // out of range
    self                          // self assignment
    } type_err;

//typedef enum {no, yes} no_yes;

