// Created by Sergey Litvinov on 08.03.2021
// Copyright 2021 ETH Zurich

#include "system.h"
#ifdef _WIN32
#include "system_windows.ipp"
#else
#include "system_unix.ipp"
#endif
