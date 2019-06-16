/**
 * @file BaseSettings.h
 */
#ifndef BASESETTINGS_H
#define BASESETTINGS_H

#include <vector>
#include <iostream>

#include "json.h"

class BaseSettings {
	public:
		BaseSettings(
			const nlohmann::json &parent_node,
			std::string section_name);

		std::vector<std::string> parse_list(
			std::string list, char delimiter) const;

	protected:
		template <typename T>
		T parse(std::string var) {
			
			T v;
			try {
				v = json_node.at(var).get<T>();
			}
			catch(std::exception &err) {
				throw std::string(
					"Error when parsing variable \"" + var + "\" in section \"" +
					section_name + "\":\n" + err.what());
			}
			return v;
		}

		template <typename T>
		std::vector<T> parse_vector(std::string var, unsigned int expected_size=0) {

			std::vector<T> v;
			try {
				v = json_node.at(var).get<std::vector<T> >();
			}
			catch(std::exception &err) {
				throw std::string(
					"Error when parsing variable \"" + var + "\" in section \"" +
					section_name + "\":\n" + err.what());
			}
			if(expected_size>0 && v.size()!=expected_size)
				throw std::string(
					"Error when parsing variable \"" + var + "\" in section \"" +
					section_name + "\": wrong array size");
			return v;
		}

		nlohmann::json json_node;
		std::string section_name;
};

#endif
