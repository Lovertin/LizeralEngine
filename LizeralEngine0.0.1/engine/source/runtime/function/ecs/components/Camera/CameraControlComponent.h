#pragma once
#include "runtime/function/ecs/components/component.h"
namespace Lizeral{

    enum class CameraMoveType {
        FreeFly,    // W 键向"面朝向"移动（可以飞天遁地）
        Ground      // W 键向"水平面朝向"移动（不能飞，高度靠跳跃）
    };

    REFLECTION_TYPE(CameraControlComponent)
    CLASS(CameraControlComponent:public Component,WhiteListFields){
        REFLECTION_BODY(CameraControlComponent)
    public:
        //setter
        void setMoveSpeed(float speed){ move_speed=speed; }
        void setSpeedMutiplier(float shiftspeed) { m_speedMultiplier=shiftspeed; }
        void setSensitivityX(float senX) { m_sensitivityX=senX; }
        void setSensitivityY(float senY) { m_sensitivityY=senY; }
        void setMoveType(CameraMoveType type) { m_moveType=type; }
        //这里y轴上下限在90°以内
        void setPitch(float pitch) {
            if (pitch > 89.0f) m_pitch = 89.0f;
            else if (pitch < -89.0f) m_pitch = -89.0f;
            else m_pitch = pitch;
        }
        void setYaw(float yaw) { m_yaw=yaw; }

        //getter
        float getMoveSpeed() const {return move_speed;}
        float getSpeedMutipier() const {return m_speedMultiplier;}
        float getSensitivityX() const {return m_sensitivityX;}
        float getSensitivityY() const {return m_sensitivityY;}
        CameraMoveType getMoveType() const {return m_moveType;}
        float getPitch() const {return m_pitch;}
        float getYaw() const {return m_yaw;}

        //叠加角度
        void addYaw(float delta) { setYaw(m_yaw + delta); }
        void addPitch(float delta) { setPitch(m_pitch + delta); }

    private:
        META(Enable)
        float move_speed{1.0f};

        META(Enable)
        float m_speedMultiplier{2.0f}; //按住shift加速以后的速度

        META(Enable)
        float m_sensitivityX{0.1f}; //rotate 水平灵敏度

        META(Enable)
        float m_sensitivityY{0.1f}; //rotate 垂直灵敏度

        CameraMoveType m_moveType{CameraMoveType::FreeFly}; //默认自由模式

        META(Enable)
        float m_pitch { 0.0f };

        META(Enable)
        float m_yaw { -90.0f };
    };
}