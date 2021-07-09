#pragma once

#include "SkunkBase.hh"

#include <glm/glm.hpp>

/**
 * Edge (half-open interval) in 2D space.
 */
struct mhd_e2d final
{
    using mhd_vec2_t = skunk::vec2_t;
    mhd_vec2_t p1{};
    mhd_vec2_t p2{};

public:
    mhd_vec2_t vec() const
    {
        return p2 - p1;
    }

public:
    static real_t dot(const mhd_e2d& e1, const mhd_e2d& e2)
    {
        mhd_vec2_t p1{e1.vec()};
        mhd_vec2_t p2{e2.vec()};
        return p1.inner_prod(p2);
    }
    static real_t len(const mhd_e2d& e1)
    {
        mhd_vec2_t p1{e1.vec()};
        return (p1).len();
    }

    static real_t dist(const mhd_e2d& e, const mhd_vec2_t& p)
    {
        mhd_e2d e1{p, e.p1};
        mhd_e2d e2{p, e.p2};
        return det(e2, e1) / len(e);
    }

public:
    static real_t det(const mhd_e2d& e1, const mhd_e2d& e2)
    {
        mhd_vec2_t p1{e1.vec()};
        mhd_vec2_t p2{e2.vec()};
        return p1.det(p2);
    }

public:
    static mhd_vec2_t normal(const mhd_e2d& e)
    {
        mhd_vec2_t p = e.vec();
        mhd_vec2_t n{p.y, -p.x};
        real_t l{n.len()};
        if (l != 0.0) {
            return n / l;
        } else {
            std::cerr << "Warning: normal to singular edge." << std::endl;
            return {1.0, 0.0};
        }
    }
    static skunk::vec3_t normal3(const mhd_e2d& e) {
        const auto n{ normal(e) };
        return { n.x, n.y, 0.0 };
    }

};	// struct mhd_e2d

typedef mhd_e2d mhd_seg2_t;

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Edge (half-open interval) in 2D space.
 */
struct mhd_e3d final
{
    using mhd_vec3_t = skunk::vec3_t;
    mhd_vec3_t p1{};
    mhd_vec3_t p2{};

public:
    
    mhd_vec3_t vec() const
    {
        return p2 - p1;
    }

public:
    
    static real_t dot(const mhd_e3d& e1, const mhd_e3d& e2)
    {
        mhd_vec3_t p1{e1.vec()};
        mhd_vec3_t p2{e2.vec()};
        return p1.inner_prod(p2);
    }
    
    static real_t len(const mhd_e3d& e1)
    {
        mhd_vec3_t p1{e1.vec()};
        return (p1).len();
    }

public:
    
    static mhd_vec3_t cross(const mhd_e3d& e1, const mhd_e3d& e2)
    {
        mhd_vec3_t p1{e1.vec()};
        mhd_vec3_t p2{e2.vec()};
        return p1.cross(p2);
    }

public:
    
    static mhd_vec3_t normal(const mhd_e3d& e1, const mhd_e3d& e2)
    {
        mhd_vec3_t n{cross(e1, e2)};
        real_t l{(n).len()};
        if (l != 0.0) {
            return n / l;
        } else {
            std::cerr << "Warning: normal to singular edges." << std::endl;
            return {1.0, 0.0};
        }
    }

public:
    static std::string str(const mhd_e3d& e);
};  // struct mhd_e3d
