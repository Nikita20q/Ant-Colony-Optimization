#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <limits>
#include <chrono>
#include <iomanip>
#include <SFML/Graphics.hpp>

// --- Вспомогательные классы Graph и Matrix ---
template <typename T>
class Graph
{
public:
    Graph(size_t num_vertices) : adj_matrix(num_vertices, std::vector<T>(num_vertices, 0)), num_vertices_(num_vertices) {}

    void add_edge(size_t from, size_t to, T weight)
    {
        adj_matrix[from][to] = weight;
        adj_matrix[to][from] = weight; // Предполагаем, что граф неориентированный
    }

    T get_weight(size_t from, size_t to) const
    {
        return adj_matrix[from][to];
    }

    size_t get_num_vertices() const
    {
        return num_vertices_;
    }

    const std::vector<std::vector<T>> &get_adj_matrix() const
    {
        return adj_matrix;
    }

    void print_graph() const
    {
        std::cout << "Graph Adjacency Matrix:\n";

        std::cout << std::setw(6) << " ";
        for (size_t j = 0; j < num_vertices_; ++j)
        {
            std::cout << std::setw(6) << j << " ";
        }
        std::cout << "\n";

        for (size_t i = 0; i < num_vertices_; ++i)
        {
            std::cout << std::setw(6) << i << " ";

            for (size_t j = 0; j < num_vertices_; ++j)
            {
                double weight = adj_matrix[i][j];
                if (weight > 999.0)
                {
                    std::cout << std::setw(6) << "--" << " ";
                }
                else
                {
                    std::cout << std::fixed << std::setprecision(2) << std::setw(6) << weight << " ";
                }
            }
            std::cout << "\n";
        }
    }

private:
    std::vector<std::vector<T>> adj_matrix;
    size_t num_vertices_;
};

template <typename T>
class Matrix
{
public:
    Matrix(size_t rows, size_t cols, T initial_value = 0.0) : data(rows, std::vector<T>(cols, initial_value)), rows_(rows), cols_(cols) {}

    T get(size_t row, size_t col) const
    {
        return data[row][col];
    }

    void set(size_t row, size_t col, T value)
    {
        data[row][col] = value;
    }

    size_t get_rows() const { return rows_; }
    size_t get_cols() const { return cols_; }

    T &operator()(size_t row, size_t col)
    {
        return data[row][col];
    }

    const T &operator()(size_t row, size_t col) const
    {
        return data[row][col];
    }

private:
    std::vector<std::vector<T>> data;
    size_t rows_;
    size_t cols_;
};

// --- Структуры и классы ACO ---

struct AntPath
{
    std::vector<std::size_t> vertices;
    double distance = 0;
};

struct Ant
{
    explicit Ant(std::size_t start_vertex = 0) : start_location(start_vertex), current_location(start_vertex)
    {
        path.vertices.push_back(start_vertex);
        visited.push_back(start_vertex);
    }

    AntPath path;
    std::vector<std::size_t> visited;
    std::size_t start_location, current_location;
    bool can_continue = true;

    void MakeChoice(const Graph<double> &g, const Matrix<double> &p, double a, double b); // alpha - важность феромонов, beta - важность длины пути
    double getRandomChoice();
    std::vector<std::size_t> getNeighborVertexes(const Graph<double> &g); // Куда можем двигаться
};

class AntColonyOptimization
{
public:
    explicit AntColonyOptimization(const Graph<double> &graph) : graph_(graph),
                                                                 pheromone_(graph_.get_num_vertices(), graph_.get_num_vertices(), kPheromone0_),
                                                                 ants_(10)
    {
        CreateAnts();
    }

    AntPath SolveSalesmansProblem();

private:
    const double kAlpha_ = 1.0;       // Влияние феромонов при выборе пути.
    const double kBeta_ = 2.0;        // Влияние расстояния (эвристики) при выборе пути.
    const double kPheromone0_ = 1.0;  // Важно, чтобы начальное значение феромона было ненулевым
    const double kQ_ = 5.0;           // Количество феромона, которое муравей оставляет на пройденном пути.
    const double kEvaporation_ = 0.2; // Коэффициент испарения феромонов.

    Graph<double> graph_;
    Matrix<double> pheromone_; // матрица феромонов - феромоны на рёбрах
    std::vector<Ant> ants_;

    void CreateAnts();
};

// Методы Ant

std::vector<std::size_t> Ant::getNeighborVertexes(const Graph<double> &g)
{
    std::vector<std::size_t> neighbors;
    for (size_t i = 0; i < g.get_num_vertices(); ++i)
    {
        if (g.get_weight(current_location, i) > 0 && std::find(visited.begin(), visited.end(), i) == visited.end()) // Если i не в visited добавляем в непосещённые
        {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

double Ant::getRandomChoice() // Генерируем рандомное число для Ant::MakeChoice
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    return dis(gen);
}

void Ant::MakeChoice(const Graph<double> &g, const Matrix<double> &p, double a, double b)
{
    std::vector<std::size_t> neighbors = getNeighborVertexes(g);

    if (neighbors.empty()) // Тупик. Пытаемся вернуться в начальный город
    {
        can_continue = false;
        if (g.get_weight(current_location, start_location) > 0)
        {
            path.vertices.push_back(start_location);
            path.distance += g.get_weight(current_location, start_location);
            current_location = start_location;
            visited.push_back(start_location);
        }
        return;
    }

    std::vector<double> probabilities; // Вероятности выбрать соседний город
    double total_probability = 0.0;    // Служит для нормализации (на неё будем делить полученные вероятности)

    for (size_t neighbor : neighbors)
    {
        double pheromone_level = p.get(current_location, neighbor);
        double distance = g.get_weight(current_location, neighbor);
        double probability = pow(pheromone_level, a) * pow(1.0 / distance, b);
        probabilities.push_back(probability);
        total_probability += probability;
    }

    // Нормализация вероятностей
    for (double &prob : probabilities)
    {
        prob /= total_probability;
    }

    // Выбор следующего города на основе вероятностей
    double random_value = getRandomChoice(); // Число от 0 до 1
    double cumulative_probability = 0.0;     // Счётчик для нормализованных вероятностей
    size_t next_city_index = 0;

    for (size_t i = 0; i < neighbors.size(); ++i)
    {
        cumulative_probability += probabilities[i];
        if (random_value <= cumulative_probability)
        {
            next_city_index = i;
            break;
        }
    }

    std::size_t next_city = neighbors[next_city_index];
    path.vertices.push_back(next_city);
    path.distance += g.get_weight(current_location, next_city);
    visited.push_back(next_city);
    current_location = next_city;
}

void AntColonyOptimization::CreateAnts()
{
    ants_.clear(); // Что бы не копить муравьёв
    for (size_t i = 0; i < graph_.get_num_vertices(); ++i)
    {
        ants_.push_back(Ant(i)); // Размещаем на каждую вершину по муравью
    }
}

AntPath AntColonyOptimization::SolveSalesmansProblem() // Лучший маршрут
{
    size_t num_iterations = 100; // Количество итераций для алгоритма
    size_t num_cities = graph_.get_num_vertices();
    double best_path_length = std::numeric_limits<double>::infinity();
    AntPath best_path;

    for (size_t iteration = 0; iteration < num_iterations; ++iteration)
    {
        // 1.  Создание путей для всех муравьев
        for (Ant &ant : ants_)
        {
            ant.path.vertices.clear();
            ant.path.distance = 0;
            ant.visited.clear();
            ant.start_location = static_cast<size_t>(&ant - &ants_[0]); // Получаем индекс муравья в векторе
            ant.current_location = ant.start_location;
            ant.visited.push_back(ant.start_location);
            ant.path.vertices.push_back(ant.start_location);
            ant.can_continue = true;

            while (ant.visited.size() < num_cities && ant.can_continue)
            {
                ant.MakeChoice(graph_, pheromone_, kAlpha_, kBeta_);
            }

            // Возврат в начальный город
            if (ant.can_continue && graph_.get_weight(ant.current_location, ant.start_location) > 0)
            {
                ant.path.vertices.push_back(ant.start_location);
                ant.path.distance += graph_.get_weight(ant.current_location, ant.start_location);
                ant.visited.push_back(ant.start_location);
            }
            else
            {
                ant.can_continue = false; // Помечаем, что муравей не смог построить полный путь
            }
        }

        // 2.  Найти лучший путь на этой итерации
        AntPath iteration_best_path;
        double iteration_best_path_length = std::numeric_limits<double>::infinity();
        for (const Ant &ant : ants_)
        {
            if (ant.can_continue && ant.path.distance < iteration_best_path_length)
            {
                iteration_best_path_length = ant.path.distance;
                iteration_best_path = ant.path;
            }
        }

        if (iteration_best_path_length < best_path_length)
        {
            best_path_length = iteration_best_path_length;
            best_path = iteration_best_path;
        }

        // Испарение феромона
        for (size_t i = 0; i < graph_.get_num_vertices(); ++i)
        {
            for (size_t j = 0; j < graph_.get_num_vertices(); ++j)
            {
                pheromone_(i, j) *= (1.0 - kEvaporation_);
            }
        }

        // Муравьи откладывают феромон
        for (const Ant &ant : ants_)
        {
            if (ant.can_continue)
            {
                double pheromone_deposit = kQ_ / ant.path.distance;
                for (size_t i = 0; i < ant.path.vertices.size() - 1; ++i)
                {
                    size_t city1 = ant.path.vertices[i];
                    size_t city2 = ant.path.vertices[i + 1];
                    pheromone_(city1, city2) += pheromone_deposit;
                    pheromone_(city2, city1) = pheromone_(city1, city2);
                }
            }
        }

        std::cout << "Итерация " << iteration + 1 << ": Best path length = " << best_path_length << std::endl;
    }

    return best_path;
}

int main()
{
    size_t num_vertices = 50;
    Graph<double> graph(num_vertices);

    // Генерация случайных расстояний
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(10.0, 100.0); // Расстояния от 10 до 100

    // Генерация случайных координат для визуализации
    std::vector<std::pair<double, double>> city_coordinates(num_vertices);
    std::uniform_real_distribution<> coord_dis(0.1, 0.9); // Координаты от 0.1 до 0.9

    for (size_t i = 0; i < num_vertices; ++i)
    {
        city_coordinates[i] = {coord_dis(gen), coord_dis(gen)};
    }

    // Создаем полный граф со случайными расстояниями
    for (size_t i = 0; i < num_vertices; ++i)
    {
        for (size_t j = i + 1; j < num_vertices; ++j)
        {
            double distance = dis(gen);
            graph.add_edge(i, j, distance);
        }
    }

    // Дополнительно рандомезируем расстояние
    for (size_t i = 0; i < num_vertices; ++i)
    {
        for (size_t j = 0; j < num_vertices; ++j)
        {
            if (i != j && graph.get_weight(i, j) > 0)
            {
                if (rand() % 100 < 10)
                {                                                      // 10% шанс сильно изменить расстояние
                    graph.add_edge(i, j, dis(gen) * (rand() % 5 + 1)); // увелииваем или уменьшаем расстояние случайным образом
                }
            }
        }
    }

    // Вывод графа в консоль
    graph.print_graph();

    AntColonyOptimization aco(graph);
    AntPath best_path = aco.SolveSalesmansProblem();

    std::cout << "Best Path: ";
    for (size_t vertex : best_path.vertices)
    {
        std::cout << vertex << " ";
    }
    std::cout << std::endl;
    std::cout << "Best Path Length: " << best_path.distance << std::endl;

    // Визуализация с использованием SFML
    sf::RenderWindow window(sf::VideoMode(800, 600), "Ant Colony Optimization TSP");
    float scale_x = 700.0f; // Масштабирование координат по X
    float scale_y = 500.0f; // Масштабирование координат по Y
    float offset_x = 50.0f; // Смещение по X
    float offset_y = 50.0f; // Смещение по X

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear(sf::Color::White);

        // Рисуем ребра (все ребра графа)
        const auto &adj_matrix = graph.get_adj_matrix();
        for (size_t i = 0; i < num_vertices; ++i)
        {
            for (size_t j = i + 1; j < num_vertices; ++j)
            { // Рисуем только один раз для каждого ребра
                if (adj_matrix[i][j] > 0)
                {
                    sf::Vertex line[] = {
                        sf::Vertex(sf::Vector2f(city_coordinates[i].first * scale_x + offset_x, city_coordinates[i].second * scale_y + offset_y)),
                        sf::Vertex(sf::Vector2f(city_coordinates[j].first * scale_x + offset_x, city_coordinates[j].second * scale_y + offset_y))};
                    line[0].color = sf::Color::Black; // Цвет ребер графа
                    line[1].color = sf::Color::Black;
                    window.draw(line, 2, sf::Lines);
                }
            }
        }

        // Рисуем города
        for (size_t i = 0; i < num_vertices; ++i)
        {
            sf::CircleShape city(5.0f);
            city.setFillColor(sf::Color::Red);
            city.setPosition(city_coordinates[i].first * scale_x + offset_x,
                             city_coordinates[i].second * scale_y + offset_y);
            window.draw(city);

            // Добавляем текст с номером города
            sf::Font font;
            if (!font.loadFromFile("arial.ttf"))
            {
                std::cerr << "Ошибка загрузки шрифта arial.ttf" << std::endl;
                return -1;
            }

            sf::Text city_number;
            city_number.setFont(font);
            city_number.setString(std::to_string(i));
            city_number.setCharacterSize(12);
            city_number.setFillColor(sf::Color::Black);
            city_number.setPosition(city_coordinates[i].first * scale_x + offset_x + 10,
                                    city_coordinates[i].second * scale_y + offset_y + 10);
            window.draw(city_number);
        }

        // Рисуем лучший путь
        for (size_t i = 0; i < best_path.vertices.size() - 1; ++i)
        {
            sf::Vertex line[] = {
                sf::Vertex(sf::Vector2f(city_coordinates[best_path.vertices[i]].first * scale_x + offset_x,
                                        city_coordinates[best_path.vertices[i]].second * scale_y + offset_y)),
                sf::Vertex(sf::Vector2f(city_coordinates[best_path.vertices[i + 1]].first * scale_x + offset_x,
                                        city_coordinates[best_path.vertices[i + 1]].second * scale_y + offset_y))};
            line[0].color = sf::Color::Blue; // Цвет ребер лучшего пути
            line[1].color = sf::Color::Blue;
            window.draw(line, 2, sf::Lines);
        }

        window.display();
    }

    return 0;
}