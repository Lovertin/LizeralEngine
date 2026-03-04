#include <QHBoxLayout>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QString>
#include <QCheckBox>
#include <QLineEdit>
#include <QSlider>
#include <QMouseEvent>
#include <QCursor>


#include <sstream>
#include <iostream>

#include "DataTypeFactory.h"
#include "core/math/vector3.h"

namespace Lizeral {

    std::unordered_map<std::string, std::unique_ptr<IPropertyDrawer>> DataTypeUIFactory::s_DrawerRegistry;
    bool DataTypeUIFactory::s_Initialized = false;

    // 1. Float Drawer
    class FloatDrawer : public IPropertyDrawer {
    public:
        QWidget* DrawProperty(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent) override {
            QDoubleSpinBox* spinBox = new QDoubleSpinBox(parent);
            spinBox->setRange(-999999.0, 999999.0); 
            spinBox->setDecimals(2);                
            spinBox->setSingleStep(1.0);          

            // 【读取初始化数据】
            float* actual_value = static_cast<float*>(accessor.get(instance));
            if (actual_value) {
                spinBox->setValue(*actual_value);
                // std::cout << "[Init] " << accessor.getFieldName() << " bound to PTR: " << instance 
                //           << " Initial Value: " << *actual_value << std::endl;
            }

            // 【写入数据并实时打印反射监控】
            QObject::connect(spinBox, &QDoubleSpinBox::valueChanged, 
                [accessor, instance](double newValue) mutable {
                    std::cout << "\n>>> [Reflection Triggered: Float] <<<" << std::endl;
                    std::cout << "  Field Name : " << accessor.getFieldName() << std::endl;
                    std::cout << "  Target PTR : " << instance << std::endl;

                    // 1. 打印写入前内存里的真实值
                    float* old_val_ptr = static_cast<float*>(accessor.get(instance));
                    if (old_val_ptr) {
                        std::cout << "  Old Memory Val: " << *old_val_ptr << std::endl;
                    } else {
                        std::cout << "  [WARNING] Memory Access Failed before set!" << std::endl;
                    }

                    // 2. 执行反射写入
                    float updated_value = static_cast<float>(newValue);
                    std::cout << "  Writing Val: " << updated_value << "..." << std::endl;
                    accessor.set(instance, &updated_value);

                    // 3. 再次读取，验证是否真正写入成功
                    float* new_val_ptr = static_cast<float*>(accessor.get(instance));
                    if (new_val_ptr) {
                        std::cout << "  New Memory Val: " << *new_val_ptr << " (Success!)" << std::endl;
                    }
                    std::cout << ">>> ------------------------------- <<<\n" << std::endl;
                }
            );

            return spinBox;
        }
    };


    // 2. Vector3 Drawer
    class Vector3Drawer : public IPropertyDrawer {
    public:
        QWidget* DrawProperty(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent) override {
            QWidget* widget = new QWidget(parent);
            QHBoxLayout* layout = new QHBoxLayout(widget);
            layout->setContentsMargins(0, 0, 0, 0);

            QDoubleSpinBox* xBox = new QDoubleSpinBox(widget);
            QDoubleSpinBox* yBox = new QDoubleSpinBox(widget);
            QDoubleSpinBox* zBox = new QDoubleSpinBox(widget);

            auto configureBox = [](QDoubleSpinBox* box) {
                box->setRange(-999999.0, 999999.0);
                box->setDecimals(3);
                box->setSingleStep(0.1);
            };
            configureBox(xBox); configureBox(yBox); configureBox(zBox);

            layout->addWidget(new QLabel("X:")); layout->addWidget(xBox);
            layout->addWidget(new QLabel("Y:")); layout->addWidget(yBox);
            layout->addWidget(new QLabel("Z:")); layout->addWidget(zBox);

            // 【读取初始化数据】
            Vector3* vec_value = static_cast<Vector3*>(accessor.get(instance));
            if (vec_value) {
                xBox->setValue(vec_value->x);
                yBox->setValue(vec_value->y);
                zBox->setValue(vec_value->z);
                // std::cout << "[Init] " << accessor.getFieldName() << " bound to PTR: " << instance 
                //           << " Initial Value: (" << vec_value->x << ", " << vec_value->y << ", " << vec_value->z << ")" << std::endl;
            }

            // 【写入数据并实时打印反射监控】
            auto updateVector3 = [accessor, instance, xBox, yBox, zBox]() mutable {
                std::cout << "\n>>> [Reflection Triggered: Vector3] <<<" << std::endl;
                std::cout << "  Field Name : " << accessor.getFieldName() << std::endl;
                std::cout << "  Target PTR : " << instance << std::endl;

                // 准备新数据
                Vector3 updated_value(
                    static_cast<float>(xBox->value()),
                    static_cast<float>(yBox->value()),
                    static_cast<float>(zBox->value())
                );
                
                std::cout << "  Writing Val: (" << updated_value.x << ", " << updated_value.y << ", " << updated_value.z << ")..." << std::endl;
                
                // 执行反射写入
                accessor.set(instance, &updated_value);

                // 验证写入
                Vector3* verify_val = static_cast<Vector3*>(accessor.get(instance));
                if (verify_val) {
                    std::cout << "  New Memory Val: (" << verify_val->x << ", " << verify_val->y << ", " << verify_val->z << ") (Success!)" << std::endl;
                } else {
                    std::cout << "  [WARNING] Memory Access Failed after set!" << std::endl;
                }
                std::cout << ">>> --------------------------------- <<<\n" << std::endl;
            };

            QObject::connect(xBox, &QDoubleSpinBox::valueChanged, updateVector3);
            QObject::connect(yBox, &QDoubleSpinBox::valueChanged, updateVector3);
            QObject::connect(zBox, &QDoubleSpinBox::valueChanged, updateVector3);

            return widget;
        }
    };

    class HeadlineStringDrawer : public IPropertyDrawer {
    public:
        QWidget* DrawProperty(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent) override {
            QLineEdit* lineEdit = new QLineEdit(parent);
            lineEdit->setMaxLength(20); // 限制 20 字符

            // 【黑魔法样式】：平时透明无边框，hover和focus时才显示边框
            lineEdit->setStyleSheet(
                "QLineEdit { border: 1px solid transparent; background: transparent; font-size: 14px; font-weight: bold; padding: 2px; }"
                "QLineEdit:hover { border: 1px solid #555555; background: #3a3a3a; }"
                "QLineEdit:focus { border: 1px solid #4CAF50; background: #2b2b2b; }"
            );

            std::string* actual_value = static_cast<std::string*>(accessor.get(instance));
            if (actual_value) lineEdit->setText(QString::fromStdString(*actual_value));

            QObject::connect(lineEdit, &QLineEdit::textChanged, 
                [accessor, instance](const QString& newValue) mutable {
                    std::string updated_value = newValue.toStdString();
                    accessor.set(instance, &updated_value);

                    Entity currentEnt = EditorSelection::Get().GetSelected();
                    emit EditorSelection::Get().OnEntityDataModified(currentEnt);
                }
            );
            return lineEdit;
        }
    };

    // 2. Bool Drawer
    class BoolDrawer : public IPropertyDrawer {
    public:
        QWidget* DrawProperty(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent) override {
            QCheckBox* checkBox = new QCheckBox(parent);
            
            bool* actual_value = static_cast<bool*>(accessor.get(instance));
            if (actual_value) checkBox->setChecked(*actual_value);

            QObject::connect(checkBox, &QCheckBox::checkStateChanged, 
                [accessor, instance](Qt::CheckState state) mutable {
                    bool updated_value = (state == Qt::Checked);
                    accessor.set(instance, &updated_value);
                }
            );
            return checkBox;
        }
    };


    // 3. 自定义支持拖拽的 SpinBox (Godot Style)
    class DragDoubleSpinBox : public QDoubleSpinBox {
        QPoint m_lastPos;
        bool m_isDragging = false;
    public:
        DragDoubleSpinBox(QWidget* parent = nullptr) : QDoubleSpinBox(parent) {
            setCursor(Qt::SizeHorCursor); // 鼠标悬浮时变成左右拖拽图标
            setButtonSymbols(QAbstractSpinBox::NoButtons); // 隐藏上下小箭头，更像 Godot
        }
    protected:
        void mousePressEvent(QMouseEvent* event) override {
            if (event->button() == Qt::LeftButton) {
                m_isDragging = true;
                m_lastPos = event->pos();
                event->accept();
            } else { QDoubleSpinBox::mousePressEvent(event); }
        }
        void mouseMoveEvent(QMouseEvent* event) override {
            if (m_isDragging) {
                // 计算鼠标水平移动了多少像素，乘以步长
                double delta = (event->pos().x() - m_lastPos.x()) * singleStep();
                setValue(value() + delta);
                m_lastPos = event->pos();
                event->accept();
            } else { QDoubleSpinBox::mouseMoveEvent(event); }
        }
        void mouseReleaseEvent(QMouseEvent* event) override {
            if (event->button() == Qt::LeftButton) { m_isDragging = false; }
            QDoubleSpinBox::mouseReleaseEvent(event);
        }
    };

    // 无限拖拽 Vector3 Drawer
    class InfiniteVector3Drawer : public IPropertyDrawer {
    public:
        QWidget* DrawProperty(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent) override {
            QWidget* widget = new QWidget(parent);
            QHBoxLayout* layout = new QHBoxLayout(widget);
            layout->setContentsMargins(0, 0, 0, 0);

            DragDoubleSpinBox* xBox = new DragDoubleSpinBox(widget);
            DragDoubleSpinBox* yBox = new DragDoubleSpinBox(widget);
            DragDoubleSpinBox* zBox = new DragDoubleSpinBox(widget);

            auto configBox = [](DragDoubleSpinBox* box) {
                box->setRange(-999999.0, 999999.0);
                box->setDecimals(3);
                box->setSingleStep(0.1);
            };
            configBox(xBox); configBox(yBox); configBox(zBox);

            layout->addWidget(new QLabel("X:")); layout->addWidget(xBox);
            layout->addWidget(new QLabel("Y:")); layout->addWidget(yBox);
            layout->addWidget(new QLabel("Z:")); layout->addWidget(zBox);

            Vector3* vec = static_cast<Vector3*>(accessor.get(instance));
            if (vec) { xBox->setValue(vec->x); yBox->setValue(vec->y); zBox->setValue(vec->z); }

            auto updateData = [accessor, instance, xBox, yBox, zBox]() mutable {
                Vector3 newVec(xBox->value(), yBox->value(), zBox->value());
                accessor.set(instance, &newVec);
            };

            QObject::connect(xBox, &QDoubleSpinBox::valueChanged, updateData);
            QObject::connect(yBox, &QDoubleSpinBox::valueChanged, updateData);
            QObject::connect(zBox, &QDoubleSpinBox::valueChanged, updateData);

            return widget;
        }
    };

    // 4. 带范围解析的 Euler 角度调节器 (底层存 Quaternion)
    class RangeQuaternionDrawer : public IPropertyDrawer {
    public:
        QWidget* DrawProperty(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent) override {
            QWidget* widget = new QWidget(parent);
            QVBoxLayout* mainLayout = new QVBoxLayout(widget);
            mainLayout->setContentsMargins(0, 0, 0, 0);

            // 【核心】：动态解析 Range MetaData！告别硬编码！
            float minVal = -360.0f;
            float maxVal = 360.0f;
            if (accessor.hasMetaData("Range")) {
                std::string rangeStr = accessor.getMetaData("Range"); // 拿到 "-360~360"
                size_t tildePos = rangeStr.find('~');
                if (tildePos != std::string::npos) {
                    minVal = std::stof(rangeStr.substr(0, tildePos));
                    maxVal = std::stof(rangeStr.substr(tildePos + 1));
                }
            }

            // 封装一个创建带滑动条的单轴 UI 的 Lambda
            auto createAxisUI = [&](const QString& label, QDoubleSpinBox*& outBox, QSlider*& outSlider) {
                QHBoxLayout* hLayout = new QHBoxLayout();
                outBox = new QDoubleSpinBox(widget);
                outSlider = new QSlider(Qt::Horizontal, widget);
                
                outBox->setRange(minVal, maxVal);
                outSlider->setRange(static_cast<int>(minVal), static_cast<int>(maxVal));

                // SpinBox 和 Slider 双向绑定
                QObject::connect(outBox, &QDoubleSpinBox::valueChanged, [outSlider](double v){ outSlider->setValue(static_cast<int>(v)); });
                QObject::connect(outSlider, &QSlider::valueChanged, [outBox](int v){ outBox->setValue(static_cast<double>(v)); });

                hLayout->addWidget(new QLabel(label));
                hLayout->addWidget(outBox);
                hLayout->addWidget(outSlider);
                mainLayout->addLayout(hLayout);
            };

            QDoubleSpinBox *xBox, *yBox, *zBox;
            QSlider *xSlider, *ySlider, *zSlider;
            createAxisUI("X", xBox, xSlider);
            createAxisUI("Y", yBox, ySlider);
            createAxisUI("Z", zBox, zSlider);

            // 读取底层数据 (Quaternion -> Euler)
            // 注意：你需要确保你的 math 库有 toEulerAngles 类似的方法
            // 如果没有，这里展示一种通用的读取占位符（你需要根据你的数学库适配）
            Quaternion* quat = static_cast<Quaternion*>(accessor.get(instance));
            if (quat) {
                float pitch = quat->getPitch().valueDegrees();  // X轴
                float yaw = quat->getYaw().valueDegrees();      // Y轴
                float roll = quat->getRoll().valueDegrees();    // Z轴

                xBox->setValue(pitch);
                yBox->setValue(yaw);
                zBox->setValue(roll);
            }

            // 写入底层数据 (Euler -> Quaternion)
            auto updateData = [accessor, instance, xBox, yBox, zBox]() mutable {
                float pitch  = static_cast<float>(xBox->value());
                float yaw = static_cast<float>(yBox->value());
                float roll = static_cast<float>(zBox->value());
                
                // 将欧拉角合成四元数 (假设你的数学库支持，或者用你 setForward 里类似的轴向乘法)
                Quaternion newQuat; 
                // newQuat.fromEulerAngles(x, y, z); // 你需要替换为你数学库对应的方法
                Radian pitchRad(pitch * Math_PI / 180.0f);
                Radian yawRad(yaw * Math_PI / 180.0f);
                Radian rollRad(roll * Math_PI / 180.0f);

                Quaternion qPitch(pitchRad, Vector3::UNIT_X);  // 绕X轴
                Quaternion qYaw(yawRad, Vector3::UNIT_Y);      // 绕Y轴
                Quaternion qRoll(rollRad, Vector3::UNIT_Z);    // 绕Z轴

                newQuat = qYaw * qPitch * qRoll;
                accessor.set(instance, &newQuat);
            };

            QObject::connect(xBox, &QDoubleSpinBox::valueChanged, updateData);
            QObject::connect(yBox, &QDoubleSpinBox::valueChanged, updateData);
            QObject::connect(zBox, &QDoubleSpinBox::valueChanged, updateData);

            return widget;
        }
    };


    // ==========================================
    // 3. 核心工厂实现
    // ==========================================
    void DataTypeUIFactory::Initialize() {
        if (s_Initialized) return;

        // 注册各种类型的绘制器
        s_DrawerRegistry["float"] = std::make_unique<FloatDrawer>();
        s_DrawerRegistry["Lizeral::Vector3"] = std::make_unique<Vector3Drawer>();
        s_DrawerRegistry["UI=Slider_Infinite_Lizeral::Vector3"] = std::make_unique<InfiniteVector3Drawer>();
        s_DrawerRegistry["bool"] = std::make_unique<BoolDrawer>();
        s_DrawerRegistry["UI=Headline_std::string"] = std::make_unique<HeadlineStringDrawer>();
        s_DrawerRegistry["UI=Slider_Range_Lizeral::Quaternion"] = std::make_unique<RangeQuaternionDrawer>();
        
        // TODO: 之后在这里注册 SliderDrawer, QuaternionDrawer 等等...

        s_Initialized = true;
    }

    QWidget* DataTypeUIFactory::CreateFieldWidget(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent) {
        if (!s_Initialized) {
            Initialize();
        }

        std::string drawerKey="";

        // // 1. 优先检查 UI Tag (高优先级策略)
        if (accessor.hasMetaData("UI")) {
            if(accessor.getMetaData("UI")=="Headline"){
                drawerKey+="UI=Headline_";
            }
            if(accessor.getMetaData("UI")=="Slider"){
                drawerKey+="UI=Slider_";
                if(accessor.hasMetaData("Range")){
                    drawerKey+="Range_";
                }
                if(accessor.hasMetaData("Infinite")){
                    drawerKey+="Infinite_";
                }
            }
            if(accessor.getMetaData("UI")=="Address"){
                drawerKey+="UI=Address_";
            }
        } 
        // 2. 退化检查数据类型

        drawerKey += accessor.getFieldTypeName();


        // 3. 在注册表中查找绘制器
        auto it = s_DrawerRegistry.find(drawerKey);
        if (it != s_DrawerRegistry.end()) {
            // 找到了！调用它的 DrawProperty
            return it->second->DrawProperty(accessor, instance, parent);
        }

        // 4. 兜底处理：如果不认识这个类型，显示一段红色的警告文字
        QLabel* fallbackLabel = new QLabel(QString("Unsupported Type: ") + QString::fromStdString(drawerKey), parent);
        fallbackLabel->setStyleSheet("color: red;");
        return fallbackLabel;
    }

} // namespace Lizeral