/**********************************************************************
shared_ptr.h - Compatibility wrapper for shared_ptr support.

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 2.
***********************************************************************/

#ifndef OB_SHARED_PTR_H
#define OB_SHARED_PTR_H

#if defined(__has_include)
  #if __has_include(<memory>)
    #include <memory>
  #elif __has_include(<tr1/memory>)
    #include <tr1/memory>
namespace std
{
  using tr1::shared_ptr;
  using tr1::weak_ptr;
  using tr1::enable_shared_from_this;
  using tr1::static_pointer_cast;
  using tr1::dynamic_pointer_cast;
  using tr1::const_pointer_cast;
}
  #endif
#else
  #if __cplusplus >= 201103L || defined(_MSC_VER)
    #include <memory>
  #else
    #include <tr1/memory>
namespace std
{
  using tr1::shared_ptr;
  using tr1::weak_ptr;
  using tr1::enable_shared_from_this;
  using tr1::static_pointer_cast;
  using tr1::dynamic_pointer_cast;
  using tr1::const_pointer_cast;
}
  #endif
#endif

#endif // OB_SHARED_PTR_H
