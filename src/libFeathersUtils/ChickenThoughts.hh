/*
 *  ______  ______   ______   ______  __  __   ______   ______   ______
 * /\  ___\/\  ___\ /\  __ \ /\__  _\/\ \_\ \ /\  ___\ /\  __ \ /\  ___\
 * \ \  __\\ \  _\  \ \  __ \\/_/\ \/\ \  __ \\ \  __\ \ \  __/ \ \___  \
 *  \ \_\   \ \_____\\ \_\ \_\  \ \_\ \ \_\ \_\\ \_____\\ \_\ \_\\/\_____\
 *   \/_/    \/_____/ \/_/\/_/   \/_/  \/_/\/_/ \/_____/ \/_/ /_/ \/_____/
 *
 * Copyright (c) 2021 Oleg Butakov
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <iostream>

namespace feathers {

/** Print a cockatiel picture. */
inline void print_cockatiel() {
    std::cout << R"(Copyright (C) Chicken Thoughts
    .#.
     -%-
      -*=
       -++-
        :+:==
         .+:.==-.
           ==...-==-.
            .==....:-==--:
              .==-.......:==-----:.
                 :==-............:-====-------::..
                    .-==-:...................:::--=====-
                        .---==-:.......................:==-
                               :=*+:......................:+-
                             .==-...........................:+
                           -+-...............................:+
                        .==:..................................-=
                      .+- .....................................*==-
                    -=-  ....................-=----==:.........-- :#+
                  ==.   ...................:*:.+**=  -*:.......:=+@@@*
               :==      ...................+-.@@@@@%  +=.......:+@@@@%
             -=:        ........-=====--:..-+.@@@@@=  *-......:==++*@-
          :+=          .......:===========-:=+--=-..=+:..........:+*.
    .-----.            .......==============-.:====-:..............-+
:===-.                  .....:===============-.......:===::.::......:*
:                       .....:================:..........:----=+*===:=:
                         ....:================:............:-==-    -++
                         .....:==============-........:----:.
                          ......-===========:......:==:
                            ......:-=====-:......:+=.
                             ...................=+
                    .:-         ...............*-
                      =                      .*.
                      -:                    :*.
                      =.                   -+
                      *                   =+
                     +:                  +=
                    +:                  *:
                   +:                  #:
                 .*:                 .*.
                -+                  .*
              :+:                  :+
            :+-                   --
)"
              << std::endl;
} // print_cockatiel

} // namespace feathers
