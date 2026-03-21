#pragma once

#include "cursor/cursor.h"

class MetaInfo
{
public:
    MetaInfo(const Cursor& cursor);

    std::string getProperty(const std::string& key) const;

    const std::unordered_map<std::string,std::string>& getProperties() const{return m_properties;}

    bool getFlag(const std::string& key) const;


private:
    typedef std::pair<std::string, std::string> Property;

    std::unordered_map<std::string, std::string> m_properties;

private:
    std::vector<Property> extractProperties(const Cursor& cursor) const;
};