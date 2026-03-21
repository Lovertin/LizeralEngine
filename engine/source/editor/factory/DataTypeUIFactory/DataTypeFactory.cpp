#include <QHBoxLayout>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QString>
#include <QCheckBox>
#include <QLineEdit>
#include <QSlider>
#include <QMouseEvent>
#include <QCursor>
#include <QPushButton>
#include <QFileDialog>
#include <QElapsedTimer>
#include <QPainter>
#include <QVariantAnimation>
#include <QPainterPath>

#include <sstream>
#include <iostream>

#include "DataTypeFactory.h"
#include "core/math/vector3.h"

namespace Lizeral {

    std::unordered_map<std::string, std::unique_ptr<IPropertyDrawer>> DataTypeUIFactory::s_DrawerRegistry;
    bool DataTypeUIFactory::s_Initialized = false;


    class DragDoubleSpinBox : public QDoubleSpinBox {
        QPoint m_lastPos;
        bool m_isDragging = false;
        QElapsedTimer m_clickTimer;

        double m_startDragRealValue = 0.0;
        
    public:

        std::function<double()> onDragStart; 
        std::function<void(double)> onDragUpdate; 
        std::function<void(double, double)> onDragEnd; 

        DragDoubleSpinBox(QWidget* parent = nullptr) : QDoubleSpinBox(parent) {
            setButtonSymbols(QAbstractSpinBox::NoButtons);
            setCursor(Qt::SizeHorCursor);
            setAlignment(Qt::AlignCenter);
            lineEdit()->setAttribute(Qt::WA_TransparentForMouseEvents, true);
            
            setStyleSheet(
                "DragDoubleSpinBox { background-color: #3d3d3d; border: 1px solid #2d2d2d; border-radius: 2px; color: #eeeeee; }"
                "DragDoubleSpinBox:hover { background-color: #4d4d4d; border: 1px solid #555555; }"
                "DragDoubleSpinBox:focus { background-color: #2b2b2b; border: 1px solid #4CAF50; color: white; }"
            );

            QObject::connect(this, &QDoubleSpinBox::editingFinished, [this]() {
                lineEdit()->setAttribute(Qt::WA_TransparentForMouseEvents, true);
            });
        }

    protected:
        void mousePressEvent(QMouseEvent* event) override {
            if (event->button() == Qt::LeftButton) {
                m_lastPos = event->globalPosition().toPoint();
                m_isDragging = false;
                m_clickTimer.start();

                if (onDragStart) {
                    m_startDragRealValue = onDragStart();
                } else {
                    m_startDragRealValue = value(); 
                }
                
                event->accept();
            } else {
                QDoubleSpinBox::mousePressEvent(event);
            }
        }

        void mouseMoveEvent(QMouseEvent* event) override {
            if (event->buttons() & Qt::LeftButton) {
                if (m_clickTimer.elapsed() > 150 || (event->globalPosition().toPoint() - m_lastPos).manhattanLength() > 3) {
                    m_isDragging = true;
                    
                    int dx = event->globalPosition().toPoint().x() - m_lastPos.x();

                    double delta = dx * singleStep();
                    double newValue = m_startDragRealValue + delta;

                    setValue(newValue);

                    if (onDragUpdate) onDragUpdate(newValue);
                    
                    setCursor(Qt::BlankCursor);
                    event->accept();
                }
            } else {
                QDoubleSpinBox::mouseMoveEvent(event);
            }
        }

        void mouseReleaseEvent(QMouseEvent* event) override {
            setCursor(Qt::SizeHorCursor);
            
            if (event->button() == Qt::LeftButton) {
                if (!m_isDragging) {
                    lineEdit()->setAttribute(Qt::WA_TransparentForMouseEvents, false);
                    setFocus();
                    lineEdit()->setFocus();
                    lineEdit()->selectAll();
                } else {
                    if (onDragEnd) {
                        onDragEnd(m_startDragRealValue, value());
                    }
                }
                m_isDragging = false;
                event->accept();
            } else {
                QDoubleSpinBox::mouseReleaseEvent(event);
            }
        }
    };

    class SliderSpinBox : public DragDoubleSpinBox {
    public:
        SliderSpinBox(QWidget* parent = nullptr) : DragDoubleSpinBox(parent) {
            lineEdit()->setStyleSheet("background: transparent; color: #c8c8c8;");
        }

    protected:
        void paintEvent(QPaintEvent* event) override {
            DragDoubleSpinBox::paintEvent(event);

            double minV = minimum();
            double maxV = maximum();
            if (maxV <= minV) return;

            double ratio = (value() - minV) / (maxV - minV);
            ratio = std::clamp(ratio, 0.0, 1.0); 

            QPainter painter(this);
            painter.setPen(Qt::NoPen);
            painter.setBrush(QColor(0, 112, 224, 100)); 
            
            QRect fillRect = rect().adjusted(1, 1, -1, -1); 
            fillRect.setWidth(static_cast<int>(fillRect.width() * ratio));
            painter.drawRect(fillRect);
        }
    };

    class AnimatedToggleSwitch : public QCheckBox {
    public:
        AnimatedToggleSwitch(QWidget* parent = nullptr) : QCheckBox(parent), m_position(isChecked() ? 1.0f : 0.0f) {
            setFixedSize(36, 18);
            setCursor(Qt::PointingHandCursor);

            QObject::connect(this, &QCheckBox::toggled, this, [this](bool checked) {
                QVariantAnimation* anim = new QVariantAnimation(this);
                anim->setDuration(150); 
                anim->setStartValue(m_position);
                anim->setEndValue(checked ? 1.0f : 0.0f);
                
                QObject::connect(anim, &QVariantAnimation::valueChanged, this, [this](const QVariant& val){
                    m_position = val.toFloat();
                    update(); 
                });
                anim->start(QAbstractAnimation::DeleteWhenStopped);
            });
        }

    protected:

        void paintEvent(QPaintEvent* event) override {
            QPainter painter(this);
            painter.setRenderHint(QPainter::Antialiasing); 

            int r = 45 + m_position * (0 - 45);
            int g = 45 + m_position * (112 - 45);
            int b = 45 + m_position * (224 - 45);
            QColor bgColor(r, g, b);

            QRectF rect(0, 0, width(), height());
            painter.setBrush(bgColor);
            painter.setPen(Qt::NoPen);
            painter.drawRoundedRect(rect, height() / 2.0, height() / 2.0);

            float circleSize = height() - 4.0f;

            float xPos = 2.0f + m_position * (width() - circleSize - 4.0f);
            
            QRectF circleRect(xPos, 2.0f, circleSize, circleSize);
            painter.setBrush(Qt::white);

            painter.setPen(QPen(QColor(0, 0, 0, 50), 1));
            painter.drawEllipse(circleRect);
        }

    private:
        float m_position;

        void checkStateSet() override {
            QCheckBox::checkStateSet();
            m_position = isChecked() ? 1.0f : 0.0f;
        }
    };


    class FloatDrawer : public IPropertyDrawer {
    public:
        QWidget* DrawProperty(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent) override {
            DragDoubleSpinBox* spinBox = new DragDoubleSpinBox(parent);
            spinBox->setRange(-999999.0, 999999.0); 
            spinBox->setDecimals(2);                
            spinBox->setSingleStep(1.0);          

            float* actual_value = static_cast<float*>(accessor.get(instance));
            if (actual_value) spinBox->setValue(*actual_value);

            Entity currentEnt = EditorSelection::Get().GetSelected();
            std::string compName = accessor.getOwnerTypeMeta().getTypeName();
            std::string fieldName = accessor.getFieldName();

            spinBox->onDragStart = [accessor, currentEnt, compName]() mutable -> double {
                void* comp = EditorContext::Get().GetComponentByName(currentEnt, compName);
                if (comp) return *static_cast<float*>(accessor.get(comp));
                return 0.0;
            };

            spinBox->onDragUpdate = [accessor, currentEnt, compName](double val) mutable {
            void* comp = EditorContext::Get().GetComponentByName(currentEnt, compName);
            
            if (comp) {
                float fVal = static_cast<float>(val);
                accessor.set(comp, &fVal);
                
                PointLightComponent* pt = static_cast<PointLightComponent*>(comp);
            } 
        };

            spinBox->onDragEnd = [currentEnt, compName, fieldName](double oldVal, double newVal) {
                if (oldVal == newVal) return; 

                auto getter = [currentEnt, compName]() -> void* {
                    return EditorContext::Get().GetComponentByName(currentEnt, compName);
                };

                EditorContext::Get().GetCommandManager()->ExecuteCommand(
                    std::make_unique<EditPropertyCommand<float>>(
                        currentEnt, compName, fieldName, 
                        static_cast<float>(oldVal), static_cast<float>(newVal), getter
                    )
                );
            };

            QObject::connect(spinBox, &QDoubleSpinBox::editingFinished, 
                [spinBox, accessor, currentEnt, compName, fieldName]() mutable {
                    void* comp = EditorContext::Get().GetComponentByName(currentEnt, compName);
                    if (!comp) return;
                    
                    float oldVal = *static_cast<float*>(accessor.get(comp));
                    float newVal = static_cast<float>(spinBox->value());
                    
                    if (oldVal != newVal) {
                        auto getter = [currentEnt, compName]() -> void* { 
                            return EditorContext::Get().GetComponentByName(currentEnt, compName); 
                        };
                        
                        EditorContext::Get().GetCommandManager()->ExecuteCommand(
                            std::make_unique<EditPropertyCommand<float>>(
                                currentEnt, compName, fieldName, oldVal, newVal, getter
                            )
                        );
                    }
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
            layout->setSpacing(4); 

            auto createAxisLayout = [&](DragDoubleSpinBox*& box, QString color) {
                QHBoxLayout* axisLayout = new QHBoxLayout();
                axisLayout->setSpacing(0); 
                axisLayout->setContentsMargins(0, 0, 0, 0);
                
                QLabel* colorBar = new QLabel();
                colorBar->setFixedSize(3, 18);
                colorBar->setStyleSheet("background-color: " + color + "; border-top-left-radius: 2px; border-bottom-left-radius: 2px;");
                
                box = new DragDoubleSpinBox(widget);
                box->setRange(-999999.0, 999999.0);
                box->setDecimals(3);
                box->setSingleStep(0.1);
                box->setMinimumWidth(40); 
                box->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
                box->setStyleSheet("QDoubleSpinBox { border-left: none; border-top-left-radius: 0px; border-bottom-left-radius: 0px; }");

                axisLayout->addWidget(colorBar);
                axisLayout->addWidget(box);
                return axisLayout;
            };

            DragDoubleSpinBox *xBox, *yBox, *zBox;
            layout->addLayout(createAxisLayout(xBox, "#c25555"));
            layout->addLayout(createAxisLayout(yBox, "#55c255"));
            layout->addLayout(createAxisLayout(zBox, "#557bc2"));

            Vector3* vec_value = static_cast<Vector3*>(accessor.get(instance));
            if (vec_value) {
                xBox->setValue(vec_value->x); yBox->setValue(vec_value->y); zBox->setValue(vec_value->z);
            }

            // write back
            auto updateVector3 = [accessor, instance, xBox, yBox, zBox]() mutable {
                Vector3 updated_value(
                    static_cast<float>(xBox->value()),
                    static_cast<float>(yBox->value()),
                    static_cast<float>(zBox->value())
                );
                accessor.set(instance, &updated_value);
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
            lineEdit->setMaxLength(20);

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

    class InfiniteVector3Drawer : public IPropertyDrawer {
    public:
        QWidget* DrawProperty(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent) override {
            QWidget* widget = new QWidget(parent);
            QHBoxLayout* layout = new QHBoxLayout(widget);
            
            // cut all Margins
            layout->setContentsMargins(0, 0, 0, 0); 
            layout->setSpacing(4); 

            QString colorX = "#c25555"; 
            QString colorY = "#55c255"; 
            QString colorZ = "#557bc2"; 

            auto createAxisLayout = [&](DragDoubleSpinBox*& box, QString color) {
                QHBoxLayout* axisLayout = new QHBoxLayout();
                axisLayout->setSpacing(0);
                axisLayout->setContentsMargins(0, 0, 0, 0);
                
                QLabel* colorBar = new QLabel();
                colorBar->setFixedSize(3, 18); 
                colorBar->setStyleSheet("background-color: " + color + "; border-top-left-radius: 2px; border-bottom-left-radius: 2px;");
                
                box = new DragDoubleSpinBox(widget);
                box->setRange(-999999.0, 999999.0);
                box->setDecimals(3);
                box->setSingleStep(0.1);
                box->setMinimumWidth(40);
                box->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
                
                box->setStyleSheet("QDoubleSpinBox { border-left: none; border-top-left-radius: 0px; border-bottom-left-radius: 0px; }");

                axisLayout->addWidget(colorBar);
                axisLayout->addWidget(box);
                return axisLayout;
            };

            DragDoubleSpinBox *xBox, *yBox, *zBox;
            layout->addLayout(createAxisLayout(xBox, colorX));
            layout->addLayout(createAxisLayout(yBox, colorY));
            layout->addLayout(createAxisLayout(zBox, colorZ));

            Vector3* vec = static_cast<Vector3*>(accessor.get(instance));
            if (vec) { xBox->setValue(vec->x); yBox->setValue(vec->y); zBox->setValue(vec->z); }

            // prepare for Command context  
            Entity currentEnt = EditorSelection::Get().GetSelected();
            std::string compName = accessor.getOwnerTypeMeta().getTypeName();
            std::string fieldName = accessor.getFieldName();

            std::shared_ptr<Vector3> dragStartVec = std::make_shared<Vector3>();

            // Wrap the single-axis binding logic in a lambda
            auto bindAxis = [&](DragDoubleSpinBox* box, int axisIndex) {
                box->onDragStart = [accessor, currentEnt, compName, axisIndex, dragStartVec]() mutable -> double {
                    void* comp = EditorContext::Get().GetComponentByName(currentEnt, compName);
                    if (comp) {
                        *dragStartVec = *static_cast<Vector3*>(accessor.get(comp));
                        if (axisIndex == 0) return dragStartVec->x;
                        if (axisIndex == 1) return dragStartVec->y;
                        return dragStartVec->z;
                    }
                    return 0.0;
                };

                box->onDragUpdate = [accessor, currentEnt, compName, axisIndex](double val) mutable {
                    void* comp = EditorContext::Get().GetComponentByName(currentEnt, compName);
                    if (comp) {
                        Vector3 newVec = *static_cast<Vector3*>(accessor.get(comp));
                        if (axisIndex == 0) newVec.x = static_cast<float>(val);
                        else if (axisIndex == 1) newVec.y = static_cast<float>(val);
                        else newVec.z = static_cast<float>(val);
                        accessor.set(comp, &newVec);
                    }
                };

                box->onDragEnd = [accessor, currentEnt, compName, fieldName, dragStartVec](double oldVal, double newVal) mutable {
                    if (oldVal == newVal) return;

                    void* comp = EditorContext::Get().GetComponentByName(currentEnt, compName);
                    if (!comp) return;
                    Vector3 finalVec = *static_cast<Vector3*>(accessor.get(comp));

                    auto getter = [currentEnt, compName]() -> void* {
                        return EditorContext::Get().GetComponentByName(currentEnt, compName);
                    };

                    EditorContext::Get().GetCommandManager()->ExecuteCommand(
                        std::make_unique<EditPropertyCommand<Vector3>>(
                            currentEnt, compName, fieldName, 
                            *dragStartVec, finalVec, getter
                        )
                    );
                };

                QObject::connect(box, &QDoubleSpinBox::editingFinished, 
                    [box, accessor, currentEnt, compName, fieldName, axisIndex]() mutable {
                        void* comp = EditorContext::Get().GetComponentByName(currentEnt, compName);
                        if (!comp) return;

                        Vector3 oldVec = *static_cast<Vector3*>(accessor.get(comp));
                        Vector3 newVec = oldVec;

                        float inputVal = static_cast<float>(box->value());
                        if (axisIndex == 0) newVec.x = inputVal;
                        else if (axisIndex == 1) newVec.y = inputVal;
                        else newVec.z = inputVal;

                        if (oldVec != newVec) { 
                            auto getter = [currentEnt, compName]() -> void* { return EditorContext::Get().GetComponentByName(currentEnt, compName); };
                            EditorContext::Get().GetCommandManager()->ExecuteCommand(
                                std::make_unique<EditPropertyCommand<Vector3>>(currentEnt, compName, fieldName, oldVec, newVec, getter)
                            );
                        }
                    }
                );
            };

            bindAxis(xBox, 0);
            bindAxis(yBox, 1);
            bindAxis(zBox, 2);

            return widget;
        }
    };

    class RangeQuaternionDrawer : public IPropertyDrawer {
    public:
        QWidget* DrawProperty(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent) override {
            QWidget* widget = new QWidget(parent);
            QVBoxLayout* mainLayout = new QVBoxLayout(widget);
            mainLayout->setContentsMargins(0, 0, 0, 0);
            mainLayout->setSpacing(4); 

            // 动态解析 Range MetaData
            float minVal = -360.0f;
            float maxVal = 360.0f;
            if (accessor.hasMetaData("Range")) {
                std::string rangeStr = accessor.getMetaData("Range"); 
                size_t tildePos = rangeStr.find('~');
                if (tildePos != std::string::npos) {
                    minVal = std::stof(rangeStr.substr(0, tildePos));
                    maxVal = std::stof(rangeStr.substr(tildePos + 1));
                }
            }

            auto createAxisUI = [&](SliderSpinBox*& outBox, QString color) {
                QHBoxLayout* hLayout = new QHBoxLayout();
                hLayout->setContentsMargins(0, 0, 0, 0);
                hLayout->setSpacing(0);

                QLabel* colorBar = new QLabel();
                colorBar->setFixedSize(3, 18);
                colorBar->setStyleSheet("background-color: " + color + "; border-top-left-radius: 2px; border-bottom-left-radius: 2px;");

                outBox = new SliderSpinBox(widget);
                outBox->setRange(minVal, maxVal);
                outBox->setDecimals(1);
                outBox->setSingleStep(1.0);
                outBox->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
                outBox->setStyleSheet("QDoubleSpinBox { border-left: none; border-top-left-radius: 0px; border-bottom-left-radius: 0px; }");

                hLayout->addWidget(colorBar);
                hLayout->addWidget(outBox);
                mainLayout->addLayout(hLayout);
            };

            SliderSpinBox *xBox, *yBox, *zBox;
            createAxisUI(xBox, "#c25555"); 
            createAxisUI(yBox, "#55c255");
            createAxisUI(zBox, "#557bc2");

            Quaternion* quat = static_cast<Quaternion*>(accessor.get(instance));
            if (quat) {
                xBox->setValue(quat->getPitch().valueDegrees());
                yBox->setValue(quat->getYaw().valueDegrees());
                zBox->setValue(quat->getRoll().valueDegrees());
            }

            auto updateData = [accessor, instance, xBox, yBox, zBox]() mutable {
                float pitch = static_cast<float>(xBox->value());
                float yaw   = static_cast<float>(yBox->value());
                float roll  = static_cast<float>(zBox->value());
                
                Quaternion newQuat; 
                Radian pitchRad(pitch * Math_PI / 180.0f);
                Radian yawRad(yaw * Math_PI / 180.0f);
                Radian rollRad(roll * Math_PI / 180.0f);

                Quaternion qPitch(pitchRad, Vector3::UNIT_X); 
                Quaternion qYaw(yawRad, Vector3::UNIT_Y);     
                Quaternion qRoll(rollRad, Vector3::UNIT_Z);   

                newQuat = qYaw * qPitch * qRoll;
                accessor.set(instance, &newQuat);
            };

            QObject::connect(xBox, &QDoubleSpinBox::valueChanged, updateData);
            QObject::connect(yBox, &QDoubleSpinBox::valueChanged, updateData);
            QObject::connect(zBox, &QDoubleSpinBox::valueChanged, updateData);

            return widget;
        }
    };

    
    class AddressStringDrawer : public IPropertyDrawer {
    public:
        QWidget* DrawProperty(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent) override {
            QWidget* widget = new QWidget(parent);
            QHBoxLayout* hLayout = new QHBoxLayout(widget);

            hLayout->setContentsMargins(0, 0, 0, 0);
            hLayout->setSpacing(4);

            std::string* strPtr = static_cast<std::string*>(accessor.get(instance));
            std::string currentPath = strPtr ? *strPtr : "";

            QLineEdit* lineEdit = new QLineEdit(QString::fromStdString(currentPath), widget);

            QPushButton* browseBtn = new QPushButton("...", widget);
            browseBtn->setFixedWidth(25);

            hLayout->addWidget(lineEdit);
            hLayout->addWidget(browseBtn);

            auto updateData = [accessor, instance, lineEdit]() mutable {
                std::string newStr = lineEdit->text().toStdString();
                accessor.set(instance, &newStr);

                // Lizeral::EditorSelection::Get().NotifyDataModified();
            };

            QObject::connect(lineEdit, &QLineEdit::editingFinished, updateData);

            QObject::connect(browseBtn, &QPushButton::clicked, [lineEdit, updateData, widget]() mutable {
                QString filePath = QFileDialog::getOpenFileName(widget, "Select Asset", "", "All Files (*.*)");
                if (!filePath.isEmpty()) {
                    lineEdit->setText(filePath);
                    updateData();
                }
            });

            return widget;
        }
    };

    void DataTypeUIFactory::Initialize() {
        if (s_Initialized) return;

        // register drawers
        s_DrawerRegistry["float"] = std::make_unique<FloatDrawer>();
        s_DrawerRegistry["Lizeral::Vector3"] = std::make_unique<Vector3Drawer>();
        s_DrawerRegistry["UI=Slider_Infinite_Lizeral::Vector3"] = std::make_unique<InfiniteVector3Drawer>();
        s_DrawerRegistry["bool"] = std::make_unique<BoolDrawer>();
        s_DrawerRegistry["UI=Headline_std::string"] = std::make_unique<HeadlineStringDrawer>();
        s_DrawerRegistry["UI=Slider_Range_Lizeral::Quaternion"] = std::make_unique<RangeQuaternionDrawer>();
        s_DrawerRegistry["UI=Address_std::string"] = std::make_unique<AddressStringDrawer>();
        
        // TODO: register coming Component

        s_Initialized = true;
    }

    QWidget* DataTypeUIFactory::CreateFieldWidget(Reflection::FieldAccessor& accessor, void* instance, QWidget* parent) {
        if (!s_Initialized) {
            Initialize();
        }

        std::string drawerKey="";

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

        drawerKey += accessor.getFieldTypeName();

        auto it = s_DrawerRegistry.find(drawerKey);
        if (it != s_DrawerRegistry.end()) {
            return it->second->DrawProperty(accessor, instance, parent);
        }

        QLabel* fallbackLabel = new QLabel(QString("Unsupported Type: ") + QString::fromStdString(drawerKey), parent);
        fallbackLabel->setStyleSheet("color: red;");
        return fallbackLabel;
    }

} // namespace Lizeral