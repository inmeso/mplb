/**
 * Copyright 2019 Sina Haeri
 *
 * Authors: See AUTHORS
 *
 * Contact: [s.haeri@ed.ac.uk and/or si.haeri@gmail.com]
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice
 *    this list of conditions and the following disclaimer in the documentation
 *    and or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/**
 * @brief This is a part of MieLAM software package. assertion_failed_msg
 * throws a full error message and terminates the application. 
 *
 * @author Sina Haeri (s.haeri@ed.ac.uk)
 * @author Soroosh Haeri (soroosh.haeri@gmail.com)
 */

#ifndef BOOST_ASSERT_MSG_FAILED_HPP
#define BOOST_ASSERT_MSG_FAILED_HPP

#define BOOST_ENABLE_ASSERT_HANDLER

#include <iostream>
#include <exception>
#include <boost/assert.hpp>
#include "ops_lib_core.h"

void boost::assertion_failed_msg(char const * expr,
                                 char const * msg,
                                 char const * function,
                                 char const * file, long line)
    {
      ops_printf("\nViolating: %s\n", expr);
      ops_printf("In function: %s\n", function);
      ops_printf("In file: %s\n", file);
      ops_printf("On Line: %s\n", line);
      ops_printf("ERROR: %s\n", msg);
      throw std::exception();
    }

#endif //BOOST_ASSERT_MSG_HPP
