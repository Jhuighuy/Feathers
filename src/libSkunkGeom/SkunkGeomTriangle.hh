#pragma once

#include "SkunkGeomEdge.hh"

/**
 * Triangle in 2D space.
 */
struct mhd_triangle2d final
{
    using mhd_vec2_t = feathers::vec2_t;
    using mhd_vec3_t = feathers::vec3_t;
    mhd_vec2_t p1{};
    mhd_vec2_t p2{};
    mhd_vec2_t p3{};

    mhd_vec2_t& node(int k)
    {
        k += 1;
        switch (k) {
            case 1:
                return p1;
            case 2:
                return p2;
            case 3:
                return p3;
            default:
                std::abort();
        }
    }

    static mhd_vec3_t barycenter3(const mhd_triangle2d& t1) {
        const mhd_vec2_t bc = (t1.p1 + t1.p2 + t1.p3)/3.0;
        return {bc.x,bc.y,0.0};
    }
    mhd_vec3_t barycenter() const {
        return barycenter3(*this);
    }

    static real_t len(const mhd_triangle2d& t1)
    {
        real_t l{};
        l += glm::length(t1.p1 - t1.p2);
        l += glm::length(t1.p2 - t1.p3);
        l += glm::length(t1.p3 - t1.p1);
        return l;
    }
    
    static real_t area(const mhd_triangle2d& t1)
    {
        real_t n{glm::determinant(feathers::mat2_t(t1.p2 - t1.p1, t1.p3 - t1.p1))};
        return 0.5 * n;
    }
    real_t area() const
    {
        return area(*this);
    }
};  // struct mhd_triangle2d

typedef mhd_triangle2d mhd_tri2_t;
typedef mhd_tri2_t geom_tri2_t;

struct geom_quad2_t {
    using mhd_vec2_t = feathers::vec2_t;
    using mhd_vec3_t = feathers::vec3_t;
    mhd_vec2_t p1{};
    mhd_vec2_t p2{};
    mhd_vec2_t p3{};
    mhd_vec2_t p4{};


public:
    
    mhd_vec2_t& node(int k)
    {
        k += 1;
        switch (k) {
            case 1:
                return p1;
            case 2:
                return p2;
            case 3:
                return p3;
            case 4:
                return p4;
            default:
                std::abort();
        }
    }

    real_t area() const
    {
        return geom_tri2_t{p1, p2, p3}.area() + geom_tri2_t{p3, p4, p1}.area();
    }
    mhd_vec3_t barycenter() const
    {
        auto a1 = geom_tri2_t{p1, p2, p3}.area(), a2 = geom_tri2_t{p3, p4, p1}.area();
        return a1/(a1+a2)*geom_tri2_t{p1, p2, p3}.barycenter() +
               a2/(a1+a2)*geom_tri2_t{p3, p4, p1}.barycenter();
    }
};

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Triangle in 3D space.
 */
struct mhd_tri3d final
{
    using mhd_vec3_t = feathers::vec3_t;
    mhd_vec3_t p1{};
    mhd_vec3_t p2{};
    mhd_vec3_t p3{};

public:

    mhd_vec3_t& node(int k)
    {
        k += 1;
        switch (k) {
            case 1:
                return p1;
            case 2:
                return p2;
            case 3:
                return p3;
            default:
                std::abort();
        }
    }
    
    static mhd_vec3_t barycenter(const mhd_tri3d& t1)
    {
        return (t1.p1 + t1.p2 + t1.p3)/3.0;
    }
    
public:
    
    static real_t area(const mhd_tri3d& t1)
    {
        mhd_vec3_t n{glm::cross(t1.p2 - t1.p1, t1.p3 - t1.p1)};
        return 0.5 * glm::length(n);
    }

public:
    
    static mhd_vec3_t normal(const mhd_tri3d& t1)
    {
        mhd_vec3_t n{glm::cross(t1.p2 - t1.p1, t1.p3 - t1.p1)};
        real_t l{glm::length(n)};
        if (l > 0) {
            n /= l;
            return n;
        } else {
            std::cerr << "Warning: normal to singular triangle." << std::endl;
            return {0.0, 0.0, 0.0};
        }
    }
};  // struct mhd_tri3d

struct mhd_tetr3d_t {
    using mhd_vec3_t = feathers::vec3_t;
    mhd_vec3_t p0, p1, p2, p3;

    mhd_vec3_t& node(int k)
    {
        switch (k) {
            case 0:
                return p0;
            case 1:
                return p1;
            case 2:
                return p2;
            case 3:
                return p3;
            default:
                std::abort();
        }
    }

public:
    mhd_vec3_t barycenter() const {
        return (p0 + p1 + p2 + p3)/4.0;
    }
    real_t volume() const {
        const auto a = p0 - p3, b = p1 - p3, c = p2 - p3;
        const auto v = 1.0/6.0*glm::dot(a, glm::cross(b, c));
        //assert(v >= 0.0);
        return std::abs(v);
    }
};

struct mhd_pyra3d_t {
    using mhd_vec3_t = feathers::vec3_t;
    mhd_vec3_t p0, p1, p2, p3;
    mhd_vec3_t p4;

    
    mhd_vec3_t& node(int k)
    {
        k += 1;
        switch (k) {
            case 1:
                return p1;
            case 2:
                return p2;
            case 3:
                return p3;
            case 4:
                return p4;
            default:
                std::abort();
        }
    }

public:
    mhd_vec3_t barycenter() const {
        const mhd_tetr3d_t tetr1{p0, p1, p3, p4}, tetr2{p1, p2, p3, p4};
        const auto v1 = tetr1.volume(), v2 = tetr2.volume();
        if((v1+v2) > 0.0) {
            return v1/(v1 + v2)*tetr1.barycenter() + v2/(v1 + v2)*tetr2.barycenter();
        } else {
            return 0.5*tetr1.barycenter() + 0.5*tetr2.barycenter();
        }
    }
    real_t volume() const {
        const mhd_tetr3d_t tetr1{p0, p1, p3, p4}, tetr2{p1, p2, p3, p4};
        const auto v1 = tetr1.volume(), v2 = tetr2.volume();
        return v1 + v2;
    }
};

struct mhd_hexa3d_t {
    using mhd_vec3_t = feathers::vec3_t;
    mhd_vec3_t p0, p1, p2, p3;
    mhd_vec3_t p4, p5, p6, p7;

public:
    mhd_vec3_t barycenter() const {
        const mhd_tetr3d_t tetr1{p0, p3, p1, p4}, tetr2{p3, p2, p1, p6},
                           tetr3{p4, p5, p6, p1}, tetr4{p4, p6, p7, p3}, tetr5{p4, p3, p1, p6};
        const auto v1 = tetr1.volume(), v2 = tetr2.volume(),
                   v3 = tetr3.volume(), v4 = tetr4.volume(), v5 = tetr5.volume();
        const auto sum_v = v1 + v2 + v3 + v4 + v5;
        assert(sum_v > 0.0);
        return v1/sum_v*tetr1.barycenter() + v2/sum_v*tetr2.barycenter() +
               v3/sum_v*tetr3.barycenter() + v4/sum_v*tetr4.barycenter() + v5/sum_v*tetr5.barycenter();
    }
    real_t volume() const {
        const mhd_tetr3d_t tetr1{ p0, p3, p1, p4 }, tetr2{ p3, p2, p1, p6 },
                           tetr3{ p4, p5, p6, p1 }, tetr4{ p4, p6, p7, p3 }, tetr5{ p4, p3, p1, p6 };
        const auto v1 = tetr1.volume(), v2 = tetr2.volume(),
                   v3 = tetr3.volume(), v4 = tetr4.volume(), v5 = tetr5.volume();
        const auto sum_v = v1 + v2 + v3 + v4 + v5;
        return sum_v;
    }
};


typedef mhd_tri3d mhd_tri3_t;
