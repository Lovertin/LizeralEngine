#include "ComponentUIFactory.h"
#include "editor/factory/DataTypeUIFactory/DataTypeFactory.h" // rely on DataType factory
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

        QGroupBox* groupBox = new QGroupBox(QString::fromStdString(meta.getTypeName()), parent);
        groupBox->setStyleSheet("QGroupBox { font-weight: bold; padding-top: 15px; margin-top: 10px; }");
        
        QVBoxLayout* groupLayout = new QVBoxLayout(groupBox);

        auto fields = meta.getFields();
        
        // std::cout << "[ComponentUIFactory] Component has " << fields.size() << " fields" << std::endl;

        int fieldIndex = 0;
        for (auto& accessor : fields) {
            QHBoxLayout* rowLayout = new QHBoxLayout();
            
            QLabel* nameLabel = new QLabel(QString::fromStdString(accessor.getFieldName()), groupBox);
            nameLabel->setMinimumWidth(80);
            rowLayout->addWidget(nameLabel);

            QWidget* fieldWidget = DataTypeUIFactory::CreateFieldWidget(accessor, instance, groupBox);
            rowLayout->addWidget(fieldWidget);

            groupLayout->addLayout(rowLayout);
        }

        std::cout << "[ComponentUIFactory] Component widget created successfully: " << groupBox << std::endl;
        return groupBox;
    }
} // namespace Lizeral
