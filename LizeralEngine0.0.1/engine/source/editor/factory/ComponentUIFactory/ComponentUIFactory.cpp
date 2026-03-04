#include "ComponentUIFactory.h"
#include "editor/factory/DataTypeUIFactory/DataTypeFactory.h" // 依赖底层的数据类型工厂
#include <QGroupBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QString>
#include <iostream>

namespace Lizeral {

    QWidget* ComponentUIFactory::CreateComponentWidget(Reflection::TypeMeta& meta, void* instance, QWidget* parent) {
        if (!instance) {
            std::cout << "[ComponentUIFactory] ERROR: Instance is null!" << std::endl;
            return nullptr;
        }
        
        if (!meta.isValid()) {
            std::cout << "[ComponentUIFactory] ERROR: Meta is invalid!" << std::endl;
            return nullptr;
        }

        // 1. 创建组件的外壳：QGroupBox (自带边框和标题)
        // 使用元数据中的类名作为标题，比如 "TransformComponent"
        QGroupBox* groupBox = new QGroupBox(QString::fromStdString(meta.getTypeName()), parent);
        groupBox->setStyleSheet("QGroupBox { font-weight: bold; padding-top: 15px; margin-top: 10px; }");
        
        // 组件内部自上而下排布
        QVBoxLayout* groupLayout = new QVBoxLayout(groupBox);

        auto fields = meta.getFields();
        
        // std::cout << "[ComponentUIFactory] Component has " << fields.size() << " fields" << std::endl;

        // 3. 遍历每个字段，交给 DataTypeUIFactory 去画
        int fieldIndex = 0;
        for (auto& accessor : fields) {
            QHBoxLayout* rowLayout = new QHBoxLayout();
            
            QLabel* nameLabel = new QLabel(QString::fromStdString(accessor.getFieldName()), groupBox);
            nameLabel->setMinimumWidth(80); // 对齐
            rowLayout->addWidget(nameLabel);

            // 【调用底层工厂】
            QWidget* fieldWidget = DataTypeUIFactory::CreateFieldWidget(accessor, instance, groupBox);
            rowLayout->addWidget(fieldWidget);

            groupLayout->addLayout(rowLayout);
        }

        std::cout << "[ComponentUIFactory] Component widget created successfully: " << groupBox << std::endl;
        return groupBox;
    }
} // namespace Lizeral
