#include "maths/maths.h"

IASCII::IASCII() {

}
bool IASCII::Open(const std::string& filename, char mode) {
    switch (mode) {
        case 'w':
            _file = std::fstream(filename, std::ios::out);
            break;
        case 'r':
            _file = std::fstream(filename, std::ios::in);
            break;
        case 'a':
            _file = std::fstream(filename, std::ios::app);
            break;
    }
    return _file.is_open();
}
void IASCII::Close() {
    _file.close();
}

void IASCII::Write(const std::string& text) {
    _file.write(text.c_str(), text.length());
}
void IASCII::Flush() {
    _file.flush();
}
std::string IASCII::ReadLine() {
    std::string result;
    std::getline(_file, result);
    return result;
}