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

#pragma once
#ifndef IMAGE_HH_
#define IMAGE_HH_

#include "SkunkBase.hh"

namespace feathers {

struct sPixel {
    byte_t r = 0, g = 0, b = 0, a = 255;
}; // struct sPixel

/**
 * 2D RGBA 8-bit Image.
 */
class cImage {
private:
    uint_t m_width = 0, m_height = 0;
    sPixel* m_pixels = nullptr;

public:
    /** Unload an image. */
    ~cImage();

    /** Init an image. */
    void init(uint_t width, uint_t height);

    /** Load an image. */
    bool load(const char* path);

    /** Store an image. */
    bool store(const char* path);

    /** Image width. */
    uint_t width() const {
        return m_width;
    }
    /** Image height. */
    uint_t height() const {
        return m_height;
    }

    /** Access a pixel. */
    FEATHERS_CONST_OVERLOAD_R(
            sPixel&, const sPixel&, operator(), (uint_t x, uint_t y), {
        FEATHERS_ASSERT((x <= m_width) && (y <= m_height));
        return m_pixels[x + y*m_width];
    })
}; // class cImage

} // feathers

#endif // ifndef IMAGE_HH_
