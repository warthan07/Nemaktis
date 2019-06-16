#include <boost/regex.hpp>
#include <iostream>

#include "BaseSettings.h"

BaseSettings::BaseSettings(
		const nlohmann::json &parent_node,
		std::string section_name) :
	section_name(section_name) {

	if(section_name=="") // root section
		json_node = parent_node;
	else {
		if(parent_node.find(section_name) != parent_node.end())
			json_node = parent_node.at(section_name);
		else
			throw std::string(
				"The section \"" + section_name + "\" is missing from the setting file");
	}
}

std::vector<std::string> BaseSettings::parse_list(
		std::string list, char delimiter) const {
	

	std::stringstream ss;
	ss.str(list);

	std::vector<std::string> items;
	std::string token, item;

	std::string empty = "";
	boost::regex whitespace_re("\\s");

	while(std::getline(ss, token, delimiter)) {
		item = boost::regex_replace(token, whitespace_re, empty);
		if(item!="")
			items.push_back(item);
	}

	return items;
}
