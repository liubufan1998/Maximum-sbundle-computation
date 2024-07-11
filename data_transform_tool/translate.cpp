#include "Utility.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

bool directoryExists(const std::string &path)
{
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

bool createDirectory(const std::string &path)
{
    if (directoryExists(path))
        return true;
    int status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if (status == 0)
    {
        std::cout << "Directory \"" << path << "\" created successfully." << std::endl;
        return true;
    }
    else
    {
        std::cout << "Failed to create directory: \"" << path << "\"." << std::endl;
        return false;
    }
}

//==============================================
string input_file, file_name;
ll n, m;
vector<pii> edges;
ui *pstart, *edge_to, *deg;

void read_graph_mtx()
{
    ifstream in(input_file);
    assert(in.is_open());
    string line;
    do
    {
        getline(in, line);
    } while (line[0] == '%');

    {
        stringstream ss(line);
        if (ss >> n >> n >> m)
        {
        }
        else
        {
            cout << input_file << " failed !!!!" << endl;
            exit(1);
        }
    }
    edges.resize(m << 1);
    unordered_map<int, int> v_map;
    int id_v = 0;
    ui idx = 0;
    for (ui i = 0; i < m; i++)
    {
        ui a, b;
        getline(in, line);
        stringstream ss(line);
        if (ss >> a >> b)
        {
            if (a == b)
                continue;
            if (!v_map.count(a))
                v_map[a] = id_v++;
            if (!v_map.count(b))
                v_map[b] = id_v++;
            a = v_map[a], b = v_map[b];
            edges[idx++] = {a, b};
            edges[idx++] = {b, a};
        }
    }
    edges.resize(idx);
    n = id_v;
    unique_pii(edges, n);
    m = edges.size();
    printf("read ok, n=%lld m= %lld\n", n, m);
    fflush(stdout);
}

void read_graph_no_suffix()
{
    freopen(input_file.c_str(), "r", stdin);
    n = read(), m = read();
    edges.resize(m * 2);
    vector<int> v_map(n + 1, -1);
    int id_v = 0;
    ui idx = 0;
    ui j = 0;
    for (ll i = 0; i < m; i++)
    {
        int a = read(), b = read();

        if (v_map[a] == -1)
            v_map[a] = id_v++;
        if (v_map[b] == -1)
            v_map[b] = id_v++;

        a = v_map[a], b = v_map[b];
        assert(a < n && b < n);
        if (a == b)
            continue;
        edges[j++] = {a, b};
        edges[j++] = {b, a};
    }
    edges.resize(j);
    unique_pii(edges, id_v);
    n = id_v;
    m = edges.size();
    printf("read ok, n=%lld m= %lld\n", n, m);
    fflush(stdout);
}

string strip(string &a)
{
    string ret = "";
    bool start = 1;
    for (char ch : a)
    {
        if (start)
        {
            if (ch == ' ' || ch == '\t')
                continue;
            start = 0;
        }
        ret += ch;
    }
    return ret;
}

bool read_two_ints(string &s, ui &a, ui &b)
{
    if (s[0] == ' ' || s[0] == '\t')
    {
        s = strip(s);
    }
    if (s.size() <= 2 || s[0] < '0' || s[0] > '9')
    {
        return false;
    }
    a = 0;
    b = 0;
    bool is_a = 1;
    for (char ch : s)
    {
        if (ch >= '0' && ch <= '9')
        {
            if (is_a)
            {
                a = a * 10 + ch - '0';
            }
            else
            {
                b = b * 10 + ch - '0';
            }
        }
        else
        {
            if (is_a)
                is_a = 0;
            else
            {
                return true;
            }
        }
    }
    if (is_a)
        return false;
    return true;
}

void read_graph_edges()
{
    ifstream in(input_file);
    assert(in.is_open());
    string line;
    do
    {
        getline(in, line);
    } while (line[0] == '%' || line[0] == '#' || line[0] == '/');
    unordered_map<int, int> v_map;
    int id_v = 0;
    do
    {
        ui a, b;
        if (read_two_ints(line, a, b))
        {
            if (!v_map.count(a))
                v_map[a] = id_v++;
            if (!v_map.count(b))
                v_map[b] = id_v++;
            a = v_map[a], b = v_map[b];
            if (a == b)
                continue;
            edges.push_back({a, b});
            edges.push_back({b, a});
        }
    } while (getline(in, line));
    n = id_v;
    unique_pii(edges, n);
    m = edges.size();
    printf("read ok, n=%lld m= %lld\n", n, m);
    fflush(stdout);
}

void read_graph()
{
    string suffix = get_file_name_suffix(input_file);
    if (suffix == "mtx")
    {
        read_graph_mtx();
        return;
    }
    else if (suffix.size() == 0)
    {
        read_graph_no_suffix();
        return;
    }
    else if (suffix == "edges")
    {
        read_graph_edges();
        return;
    }
    freopen(input_file.c_str(), "r", stdin);
    n = read(), m = read();
    edges.resize(m * 2);
    ll j = 0;
    for (ll i = 0; i < m; i++)
    {
        int a = read(), b = read();
        assert(a < n && b < n);
        if (a == b)
            continue;
        edges[j++] = {a, b};
        edges[j++] = {b, a};
    }
    edges.resize(j);
    unique_pii(edges, n);
    m = edges.size();
    printf("read ok, n=%lld m= %lld\n", n, m);
    fflush(stdout);
}

void dump_Chang(ui deg[], ui ne[])
{
    string dir_name = "./" + file_name;
    if (!createDirectory(dir_name))
    {
        printf("Failed to mkdir: %s \n", dir_name.c_str());
    }
    ui size_int = sizeof(ui);
    ui v_n = n, v_m = m;
    {
        string deg_file = dir_name + "/b_degree.bin";
        FILE *out = fopen(deg_file.c_str(), "wb");
        assert(out != nullptr);

        fwrite(&size_int, sizeof(ui), 1, out);
        fwrite(&v_n, sizeof(ui), 1, out);
        fwrite(&v_m, sizeof(ui), 1, out);
        fwrite(deg, sizeof(ui), n, out);
        fflush(out);
        fclose(out);
    }
    {
        string adj_file = dir_name + "/b_adj.bin";
        FILE *out = fopen(adj_file.c_str(), "wb");
        assert(out != nullptr);
        fwrite(ne, sizeof(ui), m, out);
        fflush(out);
        fclose(out);
    }
    cout << dir_name << " Chang ok" << endl;
}

// for degeneracy gap and kpbb
void dump_bin()
{
    ui *deg = new ui[n];
    memset(deg, 0, sizeof(ui) * n);
    for (auto &h : edges)
        deg[h.x]++;
    ui size_int = sizeof(ui);
    ui v_n = n, v_m = m;
    ui *ne = new ui[edges.size()];
    for (ll i = 0; i < m; i++)
        ne[i] = edges[i].y;
    {

        string name = "./" + file_name + ".bin";
        FILE *out = fopen(name.c_str(), "wb");
        assert(out != nullptr);
        fwrite(&size_int, sizeof(ui), 1, out);
        fwrite(&v_n, sizeof(ui), 1, out);
        fwrite(&v_m, sizeof(ui), 1, out);
        fwrite(deg, sizeof(ui), n, out);
        fwrite(ne, sizeof(ui), m, out);
        fflush(out);
        fclose(out);
        cout << name << " bin ok" << endl;
    }

    delete[] deg;
    delete[] ne;
}

int main(int argc, char *argv[])
{
    input_file = string(argv[1]);
    file_name = get_file_name_without_suffix(input_file);
    cout << "File: " << file_name << endl;
    read_graph();
    dump_bin();

    return 0;
}