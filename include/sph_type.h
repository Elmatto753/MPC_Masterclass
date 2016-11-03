/** File:		sph_type.h
 ** Author:		Dongli Zhang
 ** Contact:	dongli.zhang0129@gmail.com
 **
 ** Copyright (C) Dongli Zhang 2013
 **
 ** This program is free software;  you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License as published by
 ** the Free Software Foundation; either version 2 of the License, or
 ** (at your option) any later version.
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY;  without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
 ** the GNU General Public License for more details.
 **
 ** You should have received a copy of the GNU General Public License
 ** along with this program;  if not, write to the Free Software 
 ** Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef __SPHTYPE_H__
#define __SPHTYPE_H__

typedef unsigned int uint;

//struct float3
//{
//	float x;
//	float y;
//	float z;
//};

struct int3
{
	int m_x;
	int m_y;
	int m_z;
};

struct uint3
{
	uint m_x;
	uint m_y;
	uint m_z;
};

class vec3
{
public:
  vec3() : m_x(0.0f), m_y(0.0f), m_z(0.0f) {}

  vec3(float _x, float _y, float _z) : m_x(_x), m_y(_y), m_z(_z) {}

  struct
  {
    float m_x;
    float m_y;
    float m_z;
  };

void operator+=( const vec3& _v) { m_x += _v.m_x;
                                   m_y += _v.m_y;
                                   m_y += _v.m_z; }

void operator-=( const vec3& _v) { m_x -= _v.m_x;
                                   m_y -= _v.m_y;
                                   m_y -= _v.m_z; }

void operator*=( float _v) { m_x *= _v;
                             m_y *= _v;
                             m_y *= _v; }

void operator/=( float _v) { m_x /= _v;
                             m_y /= _v;
                             m_y /= _v; }

vec3 operator/( float _v) { return vec3(m_x/_v, m_y/_v, m_z/_v); }

vec3 operator*( float _v) { return vec3(m_x*_v, m_y*_v, m_z*_v); }

vec3 operator+( const vec3 &_v) { return vec3(m_x+_v.m_x, m_y+_v.m_y, m_z+_v.m_z); }

vec3 operator-( const vec3 &_v) { return vec3(m_x-_v.m_x, m_y-_v.m_y, m_z-_v.m_z); }

vec3 operator=( const vec3 &_v) { m_x=_v.m_x, m_y=_v.m_y, m_z=_v.m_z; return *this; }

};

#endif
