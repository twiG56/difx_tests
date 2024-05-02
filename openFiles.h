//
// Created by emichel on 01/05/24.
//

#ifndef TESTS_DIFX_OPENFILES_H
#define TESTS_DIFX_OPENFILES_H

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

void openFiles(int numantennas, char **antFiles, FILE **antStream);

int readdata(int bytestoread, FILE **antStream, uint8_t **inputdata, int numStreams);


#endif //TESTS_DIFX_OPENFILES_H
