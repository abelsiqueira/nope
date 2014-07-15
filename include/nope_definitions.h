/**
 * NOPE
 * Copyright (C) 2014  Abel Soares Siqueira
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef nope_definitions_h
#define nope_definitions_h

#define UNUSED(x) (void)(x)

#include <stdbool.h>
#include "nope_interface.h"

#define STATUS_UNINITIALIZED -1
#define STATUS_OK 0
#define STATUS_READY 1
#define STATUS_PROCESSING 2
#define STATUS_PROCESSED 3
#define STATUS_ERROR 10

#define PPEPS 1e-12

#endif
