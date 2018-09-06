#ifndef __INPUT_H__
#define __INPUT_H__

void Input(const char * prompt, int & n, int _default = 0);
void Input(const char * prompt, double & d, double _default = 0.0);
void Input(const char * prompt, char & c, char _default = 'A');
void Input(const char * prompt, char * s, char * _default = "");
void Input(const char * prompt, bool & b, bool _default);

void InputBounds(const char * prompt, int & n, int  min, int max,
                 int _default = 0);
void InputBounds(const char * prompt, double & d, double min, double max,
                 double _default = 0);

extern int InputPromptWidth;

#endif
