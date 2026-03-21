#pragma once 
#include "btBulletDynamicsCommon.h"
#include "runtime/core/meta/reflection/reflection.h"
#include <vector>
#include <algorithm>
namespace Lizeral{

class PhysicsScene{
public:
    PhysicsScene():
        m_config(0),
        m_dispatcher(0),
        m_broadphase(0),
        m_solver(0),
        m_world(0)
    {
    }
    ~PhysicsScene(){
        Cleanup();
    }

    void Initialize() {
        if (m_world) Cleanup();

        m_config = new btDefaultCollisionConfiguration();
        m_dispatcher = new btCollisionDispatcher(m_config);
        m_broadphase = new btDbvtBroadphase();
        m_solver = new btSequentialImpulseConstraintSolver();
        m_world = new btDiscreteDynamicsWorld(m_dispatcher, m_broadphase, m_solver, m_config);
        
        m_world->setGravity(btVector3(0, -10, 0));
    }

    void Cleanup() {
        if (m_world) {
            // delete Constraints
            for (int i = m_world->getNumConstraints() - 1; i >= 0; i--) {
                m_world->removeConstraint(m_world->getConstraint(i));
            }

            // delete RigidBodies
            for (int i = m_world->getNumCollisionObjects() - 1; i >= 0; i--) {
                btCollisionObject* obj = m_world->getCollisionObjectArray()[i];
                btRigidBody* body = btRigidBody::upcast(obj);
                
                // delete MotionState 
                if (body && body->getMotionState()) {
                    delete body->getMotionState();
                }
                
                m_world->removeCollisionObject(obj);
                delete obj;
            }
        }

        for (int j = 0; j < m_collisionShapes.size(); j++) {
            delete m_collisionShapes[j];
        }
        m_collisionShapes.clear();

        delete m_world; m_world = nullptr;
        delete m_solver; m_solver = nullptr;
        delete m_broadphase; m_broadphase = nullptr;
        delete m_dispatcher; m_dispatcher = nullptr;
        delete m_config; m_config = nullptr;
    }

    btDiscreteDynamicsWorld* getWorld(){
        return m_world;
    }

    void trackShape(btCollisionShape* shape) {
        m_collisionShapes.push_back(shape);
    }


private:
    btDefaultCollisionConfiguration* m_config;
    btCollisionDispatcher* m_dispatcher;
    btBroadphaseInterface* m_broadphase;
    btSequentialImpulseConstraintSolver* m_solver;
    btDiscreteDynamicsWorld* m_world;

    std::vector<btCollisionShape*> m_collisionShapes;
};
 
}
