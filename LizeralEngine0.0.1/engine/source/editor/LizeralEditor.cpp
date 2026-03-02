
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <thread>
#include <atomic>
#include <memory>

// --- Qt 头文件 ---
#include <QApplication>
#include <QMainWindow>
#include <QDockWidget>
#include <QOpenGLWidget>
#include <QVBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QCloseEvent>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QHBoxLayout>
#include <QDebug>

#include "runtime/function/ecs/registry.h"
#include "runtime/function/ecs/components/RigidBody/RigidBodyComponent.h"
#include "runtime/core/meta/serializer/serializer.h"
#include "runtime/core/meta/reflection/reflection.h"

// --- 引擎侧头文件 (假定) ---
// #include "LizeralEngine.h" 


namespace Lizeral {
    namespace Reflection {
        class TypeMetaRegister {
        public:
            static void Register(); // 声明这个静态函数存在
        };
    }
}

// =====================================================================
// 1. 视口控件 (Viewport)：未来这里将用来显示你的 OpenGL 游戏画面
// =====================================================================
class EngineViewportWidget : public QOpenGLWidget {
public:
    EngineViewportWidget(QWidget* parent = nullptr) : QOpenGLWidget(parent) {}

protected:
    // Qt 的 OpenGL 初始化
    void initializeGL() override {
        // 未来在这里初始化 GLAD 或绑定引擎的渲染上下文
        glClearColor(0.2f, 0.2f, 0.2f, 1.0f); // 深灰色背景
    }

    // Qt 的每帧渲染逻辑
    void paintGL() override {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        // 未来这里会绘制引擎传过来的 FBO 纹理
    }

    // 窗口大小改变时
    void resizeGL(int w, int h) override {
        glViewport(0, 0, w, h);
    }
};

// =====================================================================
// 2. 编辑器主窗口 (Editor Main Window)
// =====================================================================
class LizeralEditorWindow : public QMainWindow {
public:
    LizeralEditorWindow() {
        setWindowTitle("Lizeral Engine - Editor v0.0.1");
        resize(1280, 720);

        setupUI();
    }

    // 暴露一个标志位，用于告诉引擎线程何时退出
    std::atomic<bool> isEditorClosed{ false };

protected:
    // 拦截窗口关闭事件，优雅地关闭引擎线程
    void closeEvent(QCloseEvent* event) override {
        std::cout << "[Editor] Closing window, notifying engine thread..." << std::endl;
        isEditorClosed = true; 
        event->accept();
    }

private:

    Lizeral::RigidBodyComponent m_testRigidBody;

    void setupUI() {
        // --- A. 中央区域：引擎渲染视口 (Viewport) ---
        QWidget* centralWidget = new QWidget(this);
        centralWidget->setStyleSheet("background-color: #2b2b2b;"); // 深灰色背景作为占位
        QVBoxLayout* centralLayout = new QVBoxLayout(centralWidget);
        QLabel* viewportLabel = new QLabel("Engine Viewport\n(未来这里将嵌入 OpenGL 渲染画面)", centralWidget);
        viewportLabel->setAlignment(Qt::AlignCenter);
        viewportLabel->setStyleSheet("color: #888888; font-size: 18px;");
        centralLayout->addWidget(viewportLabel);
        setCentralWidget(centralWidget);

        // --- B. 左侧停靠面板：Scene Outliner (场景大纲) ---
        QDockWidget* outlinerDock = new QDockWidget("Scene Outliner", this);
        outlinerDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
        QWidget* outlinerContent = new QWidget();
        QVBoxLayout* outlinerLayout = new QVBoxLayout(outlinerContent);
        outlinerLayout->addWidget(new QLabel("未来这里显示 ECS Entity 列表..."));
        outlinerLayout->addStretch();
        outlinerDock->setWidget(outlinerContent);
        addDockWidget(Qt::LeftDockWidgetArea, outlinerDock);

        // --- C. 右侧停靠面板：Inspector (基于 TypeMeta 动态反射测试) ---
        QDockWidget* inspectorDock = new QDockWidget("Inspector (Dynamic Reflection)", this);
        inspectorDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
        QWidget* inspectorContent = new QWidget();
        QVBoxLayout* inspectorLayout = new QVBoxLayout(inspectorContent);

        // ==========================================================
        // 核心魔法：使用 TypeMeta 动态解析任何注册过的类！
        // 彻底摆脱对 .gen.h 的硬编码依赖
        // ==========================================================
        
        // 1. 获取刚体组件的全局元数据
        // 注意：如果 newFromName 获取失败导致闪退，请换成 "RigidBodyComponent" 或 "Lizeral::RigidBodyComponent"
        auto meta = Lizeral::Reflection::TypeMeta::newMetaFromName("RigidBodyComponent"); 
        
        // 2. 动态获取字段元数据 (FieldAccessor)
        auto massField = meta.getFieldByName("m_mass");
        auto kinematicField = meta.getFieldByName("m_is_kinematic");

        // ==========================================================
        // 字段 1：m_mass (float) -> 对应 Qt 的 QDoubleSpinBox
        // ==========================================================
        QHBoxLayout* massLayout = new QHBoxLayout();
        
        // 动态获取名字并创建 Label
        massLayout->addWidget(new QLabel(massField.getFieldName())); 
        
        QDoubleSpinBox* massSpinBox = new QDoubleSpinBox();
        massSpinBox->setRange(0.01, 1000.0);
        massSpinBox->setSingleStep(0.1);
        
        // 反射 Get：读取初始数据
        void* massPtr = massField.get(&m_testRigidBody);
        if (massPtr) {
            massSpinBox->setValue(*static_cast<float*>(massPtr));
        }
        
        // 反射 Set：监听 UI 修改并写回 ECS
        connect(massSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this, massField](double newValue) mutable {
            float val = static_cast<float>(newValue);
            // 连这儿都不需要知道具体算子了，直接用 field.set 动态写入！
            massField.set(&m_testRigidBody, &val);
            std::cout << "[Dynamic] Mass updated to: " << m_testRigidBody.getMass() << std::endl;
        });
        
        massLayout->addWidget(massSpinBox);
        inspectorLayout->addLayout(massLayout);

        // ==========================================================
        // 字段 2：m_is_kinematic (bool) -> 对应 Qt 的 QCheckBox
        // ==========================================================
        QHBoxLayout* kinematicLayout = new QHBoxLayout();
        kinematicLayout->addWidget(new QLabel(kinematicField.getFieldName()));
        
        QCheckBox* kinematicCheckBox = new QCheckBox();
        
        // 反射 Get：读取初始数据
        void* kinematicPtr = kinematicField.get(&m_testRigidBody);
        if (kinematicPtr) {
            kinematicCheckBox->setChecked(*static_cast<bool*>(kinematicPtr));
        }
        
        // 反射 Set：复选框点击时写回数据
        connect(kinematicCheckBox, &QCheckBox::checkStateChanged, this, [this, kinematicField](Qt::CheckState state) mutable {
            bool val = (state == Qt::Checked);
            kinematicField.set(&m_testRigidBody, &val); // 动态写入！
            std::cout << "[Dynamic] IsKinematic updated to: " 
                      << (m_testRigidBody.isKinematic() ? "True" : "False") << std::endl;
        });
        
        kinematicLayout->addWidget(kinematicCheckBox);
        inspectorLayout->addLayout(kinematicLayout);

        // ==========================================================
        // 验证按钮：随时打印当前组件的真实状态
        // ==========================================================
        QPushButton* verifyBtn = new QPushButton("打印底层 Component 数据");
        verifyBtn->setStyleSheet("margin-top: 15px; padding: 5px;");
        connect(verifyBtn, &QPushButton::clicked, this, [this]() {
            std::cout << "--- Current ECS Data ---" << std::endl;
            std::cout << "Mass: " << m_testRigidBody.getMass() << std::endl;
            std::cout << "Friction: " << m_testRigidBody.getFriction() << std::endl;
            std::cout << "IsKinematic: " << m_testRigidBody.isKinematic() << std::endl;
            std::cout << "------------------------" << std::endl;
        });
        inspectorLayout->addWidget(verifyBtn);

        inspectorLayout->addStretch(); // 把控件顶到上方
        inspectorDock->setWidget(inspectorContent);
        addDockWidget(Qt::RightDockWidgetArea, inspectorDock);
    }
};

// =====================================================================
// 3. 引擎后台线程函数
// =====================================================================
void EngineTickThread(LizeralEditorWindow* editor) {
    std::cout << "[Engine Thread] Engine started successfully!" << std::endl;

    // 引擎的死循环 (Game Loop)
    while (!editor->isEditorClosed) {
        // 1. 处理来自 UI 的 Command Queue
        // 2. ECS Tick (物理、逻辑)
        // 3. 渲染画面到 FBO
        
        // 模拟一点耗时，防止空跑把 CPU 吃满
        std::this_thread::sleep_for(std::chrono::milliseconds(16)); 
    }

    std::cout << "[Engine Thread] Engine stopped cleanly!" << std::endl;
}

// =====================================================================
// 4. 程序入口
// =====================================================================

int main(int argc, char *argv[]) {
    Lizeral::Reflection::TypeMetaRegister::Register();

    // 1. 初始化 Qt 应用程序
    QApplication app(argc, argv);

    // 2. 创建并显示编辑器主窗口
    LizeralEditorWindow editorWindow;
    editorWindow.show();

    // 3. 启动引擎线程，并把编辑器的指针传进去（方便读取退出标志）
    std::thread engineThread(EngineTickThread, &editorWindow);

    // 4. 阻塞进入 Qt 事件循环
    int exitCode = app.exec();

    // 5. 等待引擎线程安全退出
    if (engineThread.joinable()) {
        engineThread.join();
    }

    return exitCode;
}