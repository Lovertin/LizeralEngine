#pragma once

#include "runtime/core/meta/reflection/reflection.h"
#include "runtime/core/math/math.h"
#include <glm/glm.hpp>
#include <cassert>

namespace Lizeral {
    REFLECTION_TYPE(Vector2)
    CLASS(Vector2, Fields) {
        REFLECTION_BODY(Vector2)
    
    public:
        glm::vec2 m_data;
    
    public:
        // 构造函数
        Vector2():m_data(0.0f, 0.0f){}
        Vector2(float x, float y):m_data(x, y){}
        Vector2(const glm::vec2& v):m_data(v){}
        explicit Vector2(float scalar):m_data(scalar,scalar){}
        explicit Vector2(const float v[2]):m_data(v[0],v[1]){}
        
        // reflection
        float getX() const { return m_data.x; }
        void setX(float value) { m_data.x = value; }
        
        float getY() const { return m_data.y; }
        void setY(float value) { m_data.y = value; }
        
        // arithmetic operations
        Vector2 operator+(const Vector2& rhs)const{
            return m_data+rhs.m_data;
        }
        Vector2 operator-(const Vector2& rhs)const{
            return m_data-rhs.m_data;
        }
        Vector2 operator*(const float scalar)const{
            return m_data*scalar;
        }
        Vector2 operator*(const Vector2& rhs)const{
            return m_data*rhs.m_data;
        }
        Vector2 operator/(const float scalar)const{
            assert(scalar!=0);
            return m_data/scalar;
        }
        Vector2 operator/(const Vector2& rhs)const{
            return m_data/rhs.m_data;
        }
        const Vector2& operator+()const{
            return *this;
        }
        Vector2 operator-()const{
            return -m_data;
        }

        // Compound assignment
        Vector2& operator+=(const Vector2& rhs) {
            m_data += rhs.m_data;
            return *this;
        }
        Vector2& operator-=(const Vector2& rhs) {
            m_data -= rhs.m_data;
            return *this;
        }
        Vector2& operator*=(float scalar) {
            m_data *= scalar;
            return *this;
        }
        Vector2& operator*=(const Vector2& rhs){
            m_data*=rhs.m_data;
            return *this;
        }
        Vector2& operator/=(float scalar) {
            assert(scalar != 0);
            m_data/=scalar;
            return *this;
        }
        Vector2& operator/=(const Vector2& rhs){
            assert(rhs.m_data.x!=0);
            assert(rhs.m_data.y!=0);
            m_data/=rhs.m_data;
            return *this;
        }

        //friend funstions
        friend Vector2 operator+(const Vector2& lhs,float scalar){return lhs.m_data+scalar;}
        friend Vector2 operator+(float scalar,const Vector2& rhs){return rhs.m_data+scalar;}
        friend Vector2 operator-(const Vector2& lhs,float scalar){return lhs.m_data-scalar;}
        friend Vector2 operator-(float scalar,const Vector2& rhs){
            return Vector2(scalar-rhs.m_data.x,scalar-rhs.m_data.y);
        }
        friend Vector2 operator*(float scalar,const Vector2& rhs){return rhs.m_data*scalar;}
        friend Vector2 operator/(float scalar,const Vector2& rhs){
            assert(scalar!=0);
            return Vector2(scalar/rhs.m_data.x,scalar/rhs.m_data.y);
        }

        // cmp
        bool operator==(const Vector2& rhs) const {
            return m_data==rhs.m_data;
        }
        
        bool operator!=(const Vector2& rhs) const {
            return m_data!=rhs.m_data;
        }

        //copy
        float operator[](size_t i)const{
            assert(i<2);
            return (i==0?m_data.x:m_data.y);
        }
        float& operator[](size_t i){
            assert(i<2);
            return (i==0?m_data.x:m_data.y);
        }

        float length() const { return glm::length(m_data); }
        float squaredLength() const { return glm::dot(m_data, m_data); }
        
        float distance(const Vector2& rhs) const {
            return glm::distance(m_data, rhs.m_data);
        }
        
        float squaredDistance(const Vector2& rhs) const {
            return glm::dot(m_data - rhs.m_data, m_data - rhs.m_data);
        }
        
        float dotProduct(const Vector2& vec) const {
            return glm::dot(m_data, vec.m_data);
        }
        
        // 归一化函数
        float normalise() {
            float len = glm::length(m_data);
            if (len > 0.0f) {
                float inv_length = 1.0f / len;
                m_data *= inv_length;
            }
            return len;
        }
        
        Vector2 midPoint(const Vector2& vec) const {
            return Vector2((m_data.x + vec.m_data.x) * 0.5f, 
                          (m_data.y + vec.m_data.y) * 0.5f);
        }
        
        bool operator<(const Vector2& rhs) const {
            return m_data.x < rhs.m_data.x && m_data.y < rhs.m_data.y;
        }
        
        bool operator>(const Vector2& rhs) const {
            return m_data.x > rhs.m_data.x && m_data.y > rhs.m_data.y;
        }
        
        void makeFloor(const Vector2& cmp) {
            if (cmp.m_data.x < m_data.x) m_data.x = cmp.m_data.x;
            if (cmp.m_data.y < m_data.y) m_data.y = cmp.m_data.y;
        }
        
        void makeCeil(const Vector2& cmp) {
            if (cmp.m_data.x > m_data.x) m_data.x = cmp.m_data.x;
            if (cmp.m_data.y > m_data.y) m_data.y = cmp.m_data.y;
        }
        
        Vector2 perpendicular() const {
            return Vector2(-m_data.y, m_data.x);
        }
        
        float crossProduct(const Vector2& rhs) const {
            return m_data.x * rhs.m_data.y - m_data.y * rhs.m_data.x;
        }
        
        bool isZeroLength() const {
            float sqlen = (m_data.x * m_data.x) + (m_data.y * m_data.y);
            return (sqlen < (Float_EPSILON * Float_EPSILON));
        }
        
        Vector2 normalisedCopy() const {
            Vector2 ret = *this;
            ret.normalise();
            return ret;
        }
        
        Vector2 reflect(const Vector2& normal) const {
            return Vector2(*this - (2 * this->dotProduct(normal) * normal));
        }
        
        bool isNaN() const {
            return Math::isNan(m_data.x) || Math::isNan(m_data.y);
        }
        
        static Vector2 lerp(const Vector2& lhs, const Vector2& rhs, float alpha) {
            return lhs + alpha * (rhs - lhs);
        }
        
        // 静态常量声明
        static const Vector2 ZERO;
        static const Vector2 UNIT_X;
        static const Vector2 UNIT_Y;
        static const Vector2 NEGATIVE_UNIT_X;
        static const Vector2 NEGATIVE_UNIT_Y;
        static const Vector2 UNIT_SCALE;

        // translate glm
        const glm::vec2& toGlm() const { return m_data; }
        glm::vec2& toGlm() { return m_data; }
    };
    
    // 静态常量定义
    inline const Vector2 Vector2::ZERO = Vector2(0.0f, 0.0f);
    inline const Vector2 Vector2::UNIT_X = Vector2(1.0f, 0.0f);
    inline const Vector2 Vector2::UNIT_Y = Vector2(0.0f, 1.0f);
    inline const Vector2 Vector2::NEGATIVE_UNIT_X = Vector2(-1.0f, 0.0f);
    inline const Vector2 Vector2::NEGATIVE_UNIT_Y = Vector2(0.0f, -1.0f);
    inline const Vector2 Vector2::UNIT_SCALE = Vector2(1.0f, 1.0f);
}
