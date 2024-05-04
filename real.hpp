#ifndef REAL_H
#define REAL_H

// some .obj files contain doubles themselves, so using float
// would corrupt the reader.
typedef double real;

#endif
