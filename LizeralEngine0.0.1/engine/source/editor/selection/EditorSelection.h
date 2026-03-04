// EditorSelection.h
#pragma once
#include <QObject>
#include "runtime/function/ecs/entity.h"

namespace Lizeral {

    class EditorSelection : public QObject {
        Q_OBJECT
    public:
        // 全局唯一访问点
        static EditorSelection& Get() {
            static EditorSelection* instance = new EditorSelection();
            return *instance;
        }

        // --- 单例模式的严谨防御：禁用拷贝和移动 ---
        EditorSelection(const EditorSelection&) = delete;
        EditorSelection& operator=(const EditorSelection&) = delete;
        EditorSelection(EditorSelection&&) = delete;
        EditorSelection& operator=(EditorSelection&&) = delete;

        void SelectEntity(Lizeral::Entity entity);

        Lizeral::Entity GetSelected() const { return m_SelectedEntity; }

    signals:
        void OnEntitySelected(Lizeral::Entity entity);

    private:
        // --- 单例模式的严谨防御：私有化构造函数 ---
        explicit EditorSelection(QObject* parent = nullptr) : QObject(parent) {}
        
        // 析构函数由局部静态变量的生命周期自动管理
        ~EditorSelection() override = default;

    private:
        Lizeral::Entity m_SelectedEntity{ Lizeral::null_entity };
    };
    
} // namespace Lizeral