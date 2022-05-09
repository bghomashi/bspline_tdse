#pragma once

#include <vector>
#include <map>
#include <string>
#include <unordered_map>
#include <memory>
#include "utility/timer.h"

#define ProfilerPush() Profile::Push(__FUNCTION__)
#define ProfilerPop() Profile::Pop(__FUNCTION__)



class Profiler {
public:
    typedef std::shared_ptr<Profiler> Ptr_t;
    
    virtual void Push(const std::string& name);       // does this work? I doubt it
    virtual void Pop(const std::string& name);

    virtual bool PrintTo(const std::string& filename);
    virtual void Print();
};




namespace Profile {
    static void Push(const std::string& name);       // does this work? I doubt it
    static void Pop(const std::string& name);

    static bool PrintTo(const std::string& filename);
    static void Print();
};
