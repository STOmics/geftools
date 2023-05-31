@echo off

set current_dir=%~dp0
echo CPU Architecture is %CPU_ARCHITECTURE% 
echo "generate vs project geftools"

cmake -G "Visual Studio 17 2022" -A x64 %current_dir%CMakeLists.txt -B %current_dir%platforms_project/Win64
cmake --build %current_dir%platforms_project/Win64 --config Release
