cmake_minimum_required(VERSION 3.10)
project(TransportProblem LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Источники
set(SOURCES
        src/main.cpp
        src/simplex/Canonical.cpp
        src/simplex/Canonical.cpp
        src/simplex/SimplexSolver.h
        src/ProblemData.cpp
        src/ProblemData.h
        # другие .cpp ...
)

add_executable(${PROJECT_NAME} ${SOURCES})

# include для своих заголовков
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src)

# Найти Eigen (header-only). Требуем минимум 3.3 (можно убрать версию, если не нужно)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Подключить Eigen к таргету. Это добавит include-путь к Eigen.
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)

# Опционально — полезные предупреждения
if (MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4 /permissive-)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
endif()