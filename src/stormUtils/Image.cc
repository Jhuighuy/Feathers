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

#include "Image.hh"

#define STB_IMAGE_IMPLEMENTATION 1
#define STB_IMAGE_STATIC 1
#include <stb/stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION 1
#define STB_IMAGE_WRITE_STATIC 1
#include <stb/stb_image_write.h>

namespace Storm {

/**
 * Unload an image.
 */
Image2D::~Image2D() {
    free(m_pixels);
} // Image2D::cImage2D2D

/**
 * Init an image.
 */
void Image2D::init(size_t width, size_t height, Pixel pixel) {
    FEATHERS_ASSERT(m_pixels == nullptr);
    m_width = width, m_height = height;
    m_pixels = static_cast<Pixel*>(
        malloc(size_t(m_width)*m_height*sizeof(*m_pixels)));
    std::fill_n(m_pixels, size_t(m_width)*m_height, pixel);
} // Image2D::init

/**
 * Load an image.
 */
bool Image2D::load(const char* path) {
    FEATHERS_ASSERT(m_pixels == nullptr);
    int channels_in_file;
    m_pixels = reinterpret_cast<Pixel*>(stbi_load(path,
                                                  reinterpret_cast<int*>(&m_width), reinterpret_cast<int*>(&m_height),
                                                  &channels_in_file, STBI_rgb_alpha));
    return m_pixels != nullptr;
} // Image2D::load

/**
 * Store an image.
 */
bool Image2D::store(const char* path) {
    const char* ext = path + strlen(path);
    while (*ext != '.' && ext != path) ext -= 1;
    int result = 1;
    if (strcmp(ext, ".png") == 0) {
        result = stbi_write_png(
            path, int(m_width), int(m_height), STBI_rgb_alpha, m_pixels, int(4*m_width));
    } else if (strcmp(ext, ".bmp") == 0) {
        result = stbi_write_bmp(
            path, int(m_width), int(m_height), STBI_rgb_alpha, m_pixels);
    } else if (strcmp(ext, ".tga") == 0) {
        result = stbi_write_tga(
            path, int(m_width), int(m_height), STBI_rgb_alpha, m_pixels);
    } else if (strcmp(ext, ".jpg") == 0) {
        result = stbi_write_jpg(
            path, int(m_width), int(m_height), STBI_rgb_alpha, m_pixels, 95);
    }
    return result == 0;
} // Image2D::store

} // feathers
