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

/**
 * RGBA
 */
struct sPixel {
    union {
        size_t rgba;
        struct {
            byte_t r, g, b, a;
        };
    };

    constexpr sPixel(size_t rgba = 0):
        rgba(rgba) {
    }
    constexpr sPixel(byte_t r, byte_t g, byte_t b, byte_t a = 255):
        r(r), g(g), b(b), a(a) {
    }
}; // struct sPixel

static constexpr bool operator<(sPixel a, sPixel b) {
    return a.rgba < b.rgba;
}

static constexpr sPixel eBlackPixel(000, 000, 000, 255);
static constexpr sPixel eWhitePixel(255, 255, 255, 255);

static constexpr sPixel   eRedPixel(255, 000, 000, 255);
static constexpr sPixel eGreenPixel(000, 255, 000, 255);
static constexpr sPixel  eBluePixel(000, 000, 255, 255);

/**
 * 2D RGBA 8-bit Image.
 */
class cImage2D {
private:
    size_t m_width = 0, m_height = 0;
    sPixel* m_pixels = nullptr;

public:
    /** Unload an image. */
    ~cImage2D();

    /** Init an image. */
    void init(size_t width, size_t height, sPixel pixel = {});

    /** Load an image. */
    bool load(const char* path);

    /** Store an image. */
    bool store(const char* path);

    /** Image width. */
    size_t width() const {
        return m_width;
    }
    /** Image height. */
    size_t height() const {
        return m_height;
    }

    /** Access a pixel. */
    FEATHERS_CONST_OVERLOAD_R(
            sPixel&, const sPixel&, operator(), (uint_t x, uint_t y), {
        FEATHERS_ASSERT((x < m_width) && (y < m_height));
        y = m_height - 1 - y;
        return m_pixels[x + y*m_width];
    })
}; // class cImage2D

} // feathers

#endif // ifndef IMAGE_HH_
