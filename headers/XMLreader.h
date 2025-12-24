#ifndef XMLREADER_H
#define XMLREADER_H

#include "common.h"

#include <cctype>
#include <string_view>


struct xmlNode {

    std::string name;
    std::string text;
    std::unordered_map<std::string, std::string> attribute;
    std::vector<xmlNode> children;

    const xmlNode* child(std::string_view n) const {
        for (auto& c : children) if (c.name == n) return& c;
        return nullptr;
    }
};

class Scanner {

    private:
    std::string data;
    size_t i = 0;

    public:
    explicit Scanner(std::string s) : data(std::move(s)) {}

    bool eof() const {return i>= data.size();}
    char peek() const {return eof() ? '\0' : data[i]; }
    char get() {return eof() ? '\0' : data[i++];}

    bool starts_with(std::string_view sv) const {
        return data.compare(i, sv.size(), sv) == 0;
    }

    void expect(char c) {
        char got = get();
        if (got != c) {
            throw std::runtime_error(std::string("Expected '") + c + "', got '" + got + "'" );
        }
    }

    // Skip spaces
    void skip_ws() {
        while (!eof() && std::isspace(static_cast<unsigned char>(peek()))) {
            i++;
        }
    }

    std::string read_name() {
        skip_ws();
        if (eof()) throw std::runtime_error("Unexpected EOF reading name");

        auto is_start = [](char ch) {
            return std::isalpha(static_cast<unsigned char>(ch)) || ch == '_';
        };

        auto is_body = [](char ch) {
            return std::isalnum(static_cast<unsigned char>(ch)) || ch == '_' || ch == '-'|| ch == '.';
        };

        if (!is_start(peek())) {
            throw std::runtime_error("Bad name start character");
        }

        std::string out;
        out.push_back(get());

        while (!eof() && is_body(peek())) {
            out.push_back(get());
        }

        return out;
    }

    std::string read_until(char delim) {
        std:: string out;
        while (!eof() && peek() != delim) {
            out.push_back(get());
        }
        return out;
    }
};

static void trim_inplace(std::string& s) {
    auto is_ws = [](unsigned char c) {return std::isspace(c); };

    while (!s.empty() && is_ws(static_cast<unsigned char>(s.front()))) {
        s.erase(s.begin());
    }

    while (!s.empty() && is_ws(static_cast<unsigned char>(s.back()))) {
        s.pop_back();
    }
}

static std::string read_quoted(Scanner& sc) {
    sc.skip_ws();

    char q = sc.get();

    if (q != '"' && q != '\'') {
        throw std::runtime_error("Expected quote from attribute value");
    }

    std::string val = sc.read_until(q);
    sc.expect(q);

    return val;
}

static xmlNode parse_element(Scanner& sc);

static void skip_comment(Scanner& sc) {
    // Consume "<!--"
    sc.expect('<');
    sc.expect('!');
    sc.expect('-');
    sc.expect('-');

    // Advance until we see the terminating sequence "-->"
    while (!sc.eof() && !sc.starts_with("-->")) {
        sc.get();
    }
    if (sc.eof()) {
        throw std::runtime_error("Unterminated comment");
    }

    // Consume "-->"
    sc.expect('-');
    sc.expect('-');
    sc.expect('>');
}

static void parse_text_into(xmlNode& node, Scanner& sc) {
    std::string t;

    while (!sc.eof() && sc.peek() != '<') {
        t.push_back(sc.get());
    }

    trim_inplace(t);

    if (!t.empty()) {
        if (!node.text.empty()) node.text.push_back(' ');
        node.text += t;
    }
}

static xmlNode parse_element(Scanner& sc) {
    xmlNode node;

    sc.skip_ws();
    sc.expect('<');

    if (sc.peek() == '?') {
        sc.get();
        while (!sc.eof() && !sc.starts_with("?>")) sc.get();
        if (sc.eof()) throw std::runtime_error("Unterminated processing insutrction");
        sc.expect('?');
        sc.expect('>');

        sc.skip_ws();
        return parse_element(sc);
    }

    if (sc.starts_with("!--")) {
        skip_comment(sc);
        sc.skip_ws();
        return parse_element(sc);
    }

    node.name = sc.read_name();

    while (true) {
        sc.skip_ws();
        char c = sc.peek();

        if (c == '/' || c == '>') break;

        std::string key = sc.read_name();
        sc.skip_ws();
        sc.expect('=');
        std::string val = read_quoted(sc);

        node.attribute.emplace(std::move(key), std::move(val));
    }

    sc.skip_ws();

    if (sc.peek() == '/') {
        sc.get();
        sc.expect('>');
        return node;
    }

    sc.expect('>');

    while (true) {
        if (sc.eof()) {
            throw std::runtime_error("Unexpected EOF inside <" + node.name + ">");
        }

        if (sc.peek() != '<') {
            parse_text_into(node, sc);
            continue;
        }

        if (sc.starts_with("<!--")) {
            skip_comment(sc);
            continue;
        }

        if (sc.starts_with("</")) {
            sc.expect('<');
            sc.expect('/');

            std::string close = sc.read_name();

            sc.skip_ws();
            sc.expect('>');

            if (close != node.name) {
                throw std::runtime_error("Mismatched close tag </" + close + "> for <" + node.name + ">");
            }

            break;
        }

        xmlNode child = parse_element(sc);
        node.children.push_back(std::move(child));
    }

    return node;
}

static xmlNode parse_document(const std::string& xml) {
    Scanner sc(xml);

    sc.skip_ws();

    xmlNode root = parse_element(sc);

    sc.skip_ws();
    return root;
}


static std::string read_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open file: " + path);

    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

static double to_double(const std::string& s) {
    char* end = nullptr;
    double v = std::strtod(s.c_str(), &end);
    if (end == s.c_str() || *end != '\0')
        throw std::runtime_error("Bad double: '" + s + "'" );
    return v;
}

static int to_int(const std::string& s) {
    char* end = nullptr;
    long v = std::strtol(s.c_str(), &end, 10);
    if (end == s.c_str() || *end != '\0')
        throw std::runtime_error("Bad int: '" + s + "'");
    return (int)v;
}

static const xmlNode& require_child(const xmlNode& n, std::string_view name) {
    if (auto* c = n.child(name)) 
        return *c;
    throw std::runtime_error("Missing child <" + std::string(name) + "> under <" + n.name + ">");
}

static std::string require_attr(const xmlNode& n, std::string_view key) {
    auto it = n.attribute.find(std::string(key));
    if (it == n.attribute.end())
        throw std::runtime_error("Missing attribute '" + std::string(key) + "' on <" + n.name + ">");
    return it->second;
}

static bool to_bool(const std::string& s) {
    if (s == "true" || s == "1") return true;
    if (s == "false" || s == "0") return false;
    throw std::runtime_error("Bad bool: '" + s + "'");
}




#endif
