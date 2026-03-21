#include "common/precompiled.h"
#include "language_types/class.h"
#include "generator/reflection_generator.h"
#include "generator/serializer_generator.h"
#include "parser.h"

#include <filesystem>
namespace fs = std::filesystem;

#define RECURSE_NAMESPACES(kind, cursor, method, namespaces) \
    { \
        if (kind == CXCursor_Namespace) \
        { \
            auto display_name = cursor.getDisplayName(); \
            if (!display_name.empty()) \
            { \
                namespaces.emplace_back(display_name); \
                method(cursor, namespaces); \
                namespaces.pop_back(); \
            } \
        } \
    }

#define TRY_ADD_LANGUAGE_TYPE(handle, container) \
    { \
        if (handle->shouldCompile()) \
        { \
            auto file = handle->getSourceFile(); \
            m_schema_modules[file].container.emplace_back(handle); \
            m_type_table[handle->m_display_name] = file; \
        } \
    }

void MetaParser::prepare(void) {}

std::string MetaParser::getIncludeFile(std::string name)
{
    auto iter = m_type_table.find(name);
    return iter == m_type_table.end() ? std::string() : iter->second;
}

MetaParser::MetaParser(const std::string project_input_file,
                       const std::string include_file_path,
                       const std::string include_path,   
                       const std::string sys_include,
                       const std::string module_name,
                       bool              is_show_errors) :
    m_project_input_file(project_input_file),
    m_source_include_file_name(include_file_path), m_index(nullptr), m_translation_unit(nullptr),
    m_sys_include(sys_include), m_module_name(module_name), m_is_show_errors(is_show_errors)
{
    m_work_paths = Utils::split(include_path, ";");

    if(!m_work_paths.empty()){
         m_generators.emplace_back(new Generator::SerializerGenerator(
            m_work_paths[0], std::bind(&MetaParser::getIncludeFile, this, std::placeholders::_1)));
         m_generators.emplace_back(new Generator::ReflectionGenerator(
            m_work_paths[0], std::bind(&MetaParser::getIncludeFile, this, std::placeholders::_1)));
    }
}

MetaParser::~MetaParser(void)
{
    for (auto item : m_generators)
    {
        delete item;
    }
    m_generators.clear();

    if (m_translation_unit)
        clang_disposeTranslationUnit(m_translation_unit);

    if (m_index)
        clang_disposeIndex(m_index);
}

void MetaParser::finish(void)
{
    for (auto generator_iter : m_generators)
    {
        generator_iter->finish();
    }
}

// 【核心修改】替换后的 parseProject 函数
bool MetaParser::parseProject()
{
    // 1. 准备输出文件
    std::fstream include_file;
    include_file.open(m_source_include_file_name, std::ios::out);
    if (!include_file.is_open())
    {
        std::cout << "Could not open the Source Include file: " << m_source_include_file_name << std::endl;
        return false;
    }

    std::cout << "Generating the Source Include file: " << m_source_include_file_name << std::endl;

    // 2. 生成 Header Guard
    std::string output_filename = Utils::getFileName(m_source_include_file_name);
    if (output_filename.empty())
    {
        output_filename = "META_INPUT_HEADER_H";
    }
    else
    {
        Utils::replace(output_filename, ".", "_");
        Utils::replace(output_filename, " ", "_");
        Utils::toUpper(output_filename);
    }
    include_file << "#ifndef __" << output_filename << "__" << std::endl;
    include_file << "#define __" << output_filename << "__" << std::endl;

    // 3. 【重点】扫描目录逻辑
    // m_project_input_file 现在是从 CMake 传过来的 runtime 目录路径
    std::string scan_root_path = m_project_input_file;
    Utils::replace(scan_root_path, '\\', '/'); 

    if (!fs::exists(scan_root_path))
    {
        std::cout << "Error: Directory does not exist: " << scan_root_path << std::endl;
        return false;
    }

    std::cout << "Scanning directory: " << scan_root_path << std::endl;

    // 递归遍历所有 .h 文件
    for (auto& dir_entry : fs::recursive_directory_iterator(scan_root_path))
    {
        if (dir_entry.is_regular_file())
        {
            auto path = dir_entry.path();
            std::string temp_string = path.string();
            Utils::replace(temp_string, '\\', '/');

            // ==========================================================
            // 【关键修复】：黑名单过滤逻辑
            // 如果文件路径中包含 "/vendor/" 或 "/3rdparty/"，直接跳过不处理！
            // 这样 libclang 就永远不会知道 Assimp 的存在，彻底杜绝命名空间污染
            // ==========================================================
            if (temp_string.find("/vendor/") != std::string::npos || 
                temp_string.find("/3rdparty/") != std::string::npos)
            {
                continue; 
            }

            if (path.extension() == ".h")
            {
                // 写入 include
                include_file << "#include  \"" << temp_string << "\"" << std::endl;
            }
        }
    }

    include_file << "#endif" << std::endl;
    include_file.close();
    return true;
}

int MetaParser::parse(void)
{
    bool parse_include_ = parseProject();
    if (!parse_include_)
    {
        std::cerr << "Parsing project file error! " << std::endl;
        return -1;
    }

    std::cerr << "Parsing the whole project..." << std::endl;
    int is_show_errors      = m_is_show_errors ? 1 : 0;
    m_index                 = clang_createIndex(true, is_show_errors);
    
    std::string pre_include = "-I";
    std::string sys_include_temp;
    if (!(m_sys_include == "*"))
    {
        sys_include_temp = pre_include + m_sys_include;
        arguments.emplace_back(sys_include_temp.c_str());
    }

    auto paths = m_work_paths;
    for (int index = 0; index < paths.size(); ++index)
    {
        paths[index] = pre_include + paths[index];
        arguments.emplace_back(paths[index].c_str());
    }

    fs::path input_path(m_source_include_file_name);
    if (!fs::exists(input_path))
    {
        std::cerr << input_path << " does not exist" << std::endl;
        return -2;
    }

    m_translation_unit = clang_createTranslationUnitFromSourceFile(
        m_index, m_source_include_file_name.c_str(), static_cast<int>(arguments.size()), arguments.data(), 0, nullptr);
    auto cursor = clang_getTranslationUnitCursor(m_translation_unit);

    Namespace temp_namespace;

    buildClassAST(cursor, temp_namespace);

    temp_namespace.clear();

    return 0;
}

void MetaParser::generateFiles(void)
{
    std::cerr << "Start generate runtime schemas(" << m_schema_modules.size() << ")..." << std::endl;
    for (auto& schema : m_schema_modules)
    {
        for (auto& generator_iter : m_generators)
        {
            generator_iter->generate(schema.first, schema.second);
        }
    }

    finish();
}

void MetaParser::buildClassAST(const Cursor& cursor, Namespace& current_namespace)
{
    for (auto& child : cursor.getChildren())
    {
        auto kind = child.getKind();

        if (child.isDefinition() && (kind == CXCursor_ClassDecl || kind == CXCursor_StructDecl))
        {
            auto class_ptr = std::make_shared<Class>(child, current_namespace);
            TRY_ADD_LANGUAGE_TYPE(class_ptr, classes);
        }
        else
        {
            RECURSE_NAMESPACES(kind, child, buildClassAST, current_namespace);
        }
    }
}