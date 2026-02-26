CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

SRC_DIR = src
TARGET = solver

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^
	@echo "Сборка успешно завершена! Запустите: ./$(TARGET)"

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(SRC_DIR)/*.o $(TARGET)
	@echo "Временные файлы и исполняемый файл удалены."
