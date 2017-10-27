/*
 * File:   common.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * FIXME: remove this file
 *
 *
 */

#pragma once

#include "graph.h"

#include <map>
#include <set>

// TODO: move these somewhere else
TMask u6c4c_cycles[6];
TMask u3_inv_pm[3];

map<pair<TMask, TMask>, set<set<TMask>>> u5cdc_from_33pp;
