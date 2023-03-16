#ifndef INTERPRETER_COMPLEX_HPP
#define INTERPRETER_COMPLEX_HPP

#include "exprtk_complex.hpp"

#include <string>
#include <map>

template <typename T>
class Interpreter
{
        typedef exprtk::symbol_table<T> symbol_table_t;
        typedef exprtk::expression<T>   expression_t;
        typedef exprtk::parser<T>       parser_t;
    private:
        T m_parameter;
        T m_changeable_constant;
        std::string m_expression_string;
        symbol_table_t m_symbol_table;
        expression_t m_expression;
        std::map<std::string, T> m_symbol_map;
    public:  
        void addVariable(std::string name, T reference)
        {
            m_symbol_map.emplace(name, reference);
            m_symbol_table.add_variable(name, m_symbol_map.at(name));
            
        }
        Interpreter()
        {
            m_symbol_table.add_constant("i", T(0,1));
            m_symbol_table.add_constant("j", T(0,1));
            m_symbol_table.add_constants();
            addVariable("c", m_changeable_constant);
            addVariable("z", m_parameter);
            
        }

        void parse(std::string expression_string)
        {
            m_expression_string = expression_string;
            m_expression.register_symbol_table(m_symbol_table);
            parser_t parser;
            parser.compile(m_expression_string, m_expression);
        }

        T getValue()
        {
            return m_expression.value();
        }

        void setValue(std::string name, T value)
        {
            m_symbol_map.at(name) = value;
        }

};
#endif /*INTERPRTER_COMPLEX*/