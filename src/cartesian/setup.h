#ifndef SETUP
#define SETUP

#include "geom.h"
#include "arb.h"
#include "options.h"

void get_domain(geom_t geometry, arb_t nu_enclosure, options_t options,
                int triangle);

#endif
