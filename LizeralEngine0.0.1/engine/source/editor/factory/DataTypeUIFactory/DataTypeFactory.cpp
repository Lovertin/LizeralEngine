#include "DataTypeFactory.h"
#include <QHBoxLayout>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QString>
#include <iostream>

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
            // QObject::connect(spinBox, &QDoubleSpinBox::valueChanged, 
            //     [accessor, instance](double newValue) mutable {
            //         std::cout << "\n>>> [Reflection Triggered: Float] <<<" << std::endl;
            //         std::cout << "  Field Name : " << accessor.getFieldName() << std::endl;
            //         std::cout << "  Target PTR : " << instance << std::endl;

            //         // 1. 打印写入前内存里的真实值
            //         float* old_val_ptr = static_cast<float*>(accessor.get(instance));
            //         if (old_val_ptr) {
            //             std::cout << "  Old Memory Val: " << *old_val_ptr << std::endl;
            //         } else {
            //             std::cout << "  [WARNING] Memory Access Failed before set!" << std::endl;
            //         }

            //         // 2. 执行反射写入
            //         float updated_value = static_cast<float>(newValue);
            //         std::cout << "  Writing Val: " << updated_value << "..." << std::endl;
            //         accessor.set(instance, &updated_value);

            //         // 3. 再次读取，验证是否真正写入成功
            //         float* new_val_ptr = static_cast<float*>(accessor.get(instance));
            //         if (new_val_ptr) {
            //             std::cout << "  New Memory Val: " << *new_val_ptr << " (Success!)" << std::endl;
            //         }
            //         std::cout << ">>> ------------------------------- <<<\n" << std::endl;
            //     }
            // );

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

    // ==========================================
    // 3. 核心工厂实现
    // ==========================================
    void DataTypeUIFactory::Initialize() {
        if (s_Initialized) return;

        // 注册各种类型的绘制器
        s_DrawerRegistry["float"] = std::make_unique<FloatDrawer>();
        s_DrawerRegistry["Lizeral::Vector3"] = std::make_unique<Vector3Drawer>();
        
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
            drawerKey = "UI_" + accessor.getMetaData("UI")+"_";
        } 
        if(accessor.hasMetaData("Range")){
            drawerKey += "Range_";
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