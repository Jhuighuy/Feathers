/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///
/// Copyright (C) 2022 Oleg Butakov
///
/// Permission is hereby granted, free of charge, to any person
/// obtaining a copy of this software and associated documentation
/// files (the "Software"), to deal in the Software without
/// restriction, including without limitation the rights  to use,
/// copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the
/// Software is furnished to do so, subject to the following
/// conditions:
///
/// The above copyright notice and this permission notice shall be
/// included in all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
/// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
/// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
/// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
/// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
/// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
/// OTHER DEALINGS IN THE SOFTWARE.
/// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ///

#pragma once

#include <stormBase.hxx>
#include <stormUtils/Index.hxx>

namespace Storm {

/// @todo Move me to the base header.
template<class>
inline constexpr bool always_false = false;

class NodeTag;
class EdgeTag;
class FaceTag;
class CellTag;
template<class>
class MarkTag;

/// @brief Node index.
using NodeIndex = Index<NodeTag>;

/// @brief Edge index.
using EdgeIndex = Index<EdgeTag>;

/// @brief Face index.
using FaceIndex = Index<FaceTag>;

/// @brief Cell index.
using CellIndex = Index<CellTag>;

/// @brief Node Mark index.
using NodeMark = Index<MarkTag<NodeTag>>;

/// @brief Edge Mark index.
using EdgeMark = Index<MarkTag<EdgeTag>>;

/// @brief Face Mark index.
using FaceMark = Index<MarkTag<FaceTag>>;

/// @brief Cell Mark index.
using CellMark = Index<MarkTag<CellTag>>;

} // namespace Storm
