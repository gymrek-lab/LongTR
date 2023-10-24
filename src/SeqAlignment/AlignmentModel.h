#ifndef ALIGNMENT_MODEL_H_
#define ALIGNMENT_MODEL_H_

#include <iostream>

const unsigned int MAX_HOMOP_LEN = 15;
const double LOG_INS_TO_INS      = -1.0; // log(e^-1)
const double LOG_INS_TO_MATCH    = -0.458675; // log(1-e^-1)
const double LOG_DEL_TO_DEL      = -1.0; // log(e^-1)
const double LOG_DEL_TO_MATCH    = -0.458675; // log(1-e^-1)

const double LOG_MATCH_TO_MATCH = -0.00005800168;
const double LOG_MATCH_TO_INS = -10.448214728;
const double LOG_MATCH_TO_DEL = -10.448214728;

void init_alignment_model();
void print_alignment_model(std::ostream& out);

#endif
