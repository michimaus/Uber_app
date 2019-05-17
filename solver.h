// Copyright 2019 SD_Homework_Team
// Duican Mihnea Ionut -- 314CA 2019

#ifndef SOLVER_H_
#define SOLVER_H_
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <queue>
#include <iostream>
#include <string>
#include <tuple>
#include <iomanip>
#include <list>
#include <ctime>
#include <sstream>
#define prime_number 503
#define large_number 20000000
#define maxval 200000000
#define minval -200000000

////////////////////////////////////
// defining the data-types I needed for my implementation

// type used for sorting by the ditance coverd by the drivers
struct dr_dist {
	std::string name;
	int value;

	// simple constructor
	dr_dist() {
		this->value = -1;
		this->name = "aaa\0";
	}

	// establishing how to have the elements compared
	bool operator < (const dr_dist &other) {
		if (this->value > other.value) {
			return true;
		} else if (this->value == other.value) {
			if (this->name < other.name)
				return true;
		}
		return false;
	}

	bool operator == (const dr_dist &other) {
		if (this->name == other.name)
			return true;
		return false;
	}

	// another constructor for specific values
	dr_dist(const int &value, const std::string &name) {
		this->name = name;
		this->value = value;
	}

	// opperator used for updating
	dr_dist &operator += (const dr_dist & other) {
		this->value = this->value + other.value;
		return *this;
	}

	// assigment opperator
	dr_dist & operator = (const int &other) {
		this->value = other;
		this->name = "aaa\0";

		return *this;
	}
	// overloading the way data is printed
	friend std::ofstream &operator <<(std::ofstream &o_d, const dr_dist &data);
};
std::ofstream &operator <<(std::ofstream &o_d, const dr_dist &data) {
	o_d << data.name << ":" << data.value;
	return o_d;
}

// type used for sorting by the rating of the drivers
struct dr_star {
	std::string name;
	float num_stars;
	int num_races;
	float rating;

	dr_star() {
		this->rating = -1;
		this->name = "aaa\0";
	}

	dr_star(const int &num_races, const float &num_stars,
		const std::string &name) {
		this->name = name;
		this->num_races = num_races;
		this->num_stars = num_stars;

		if (this->num_races == 0) {
			this->rating = 0;
		} else {
			this->rating = this->num_stars / this->num_races;
		}
	}

	bool operator == (const dr_star &other) {
		if (this->name == other.name)
			return true;
		return false;
	}

	dr_star &operator += (const dr_star & other) {
		this->num_races = this->num_races + other.num_races;
		this->num_stars = this->num_stars + other.num_stars;

		this->rating = this->num_stars / this->num_races;
		return *this;
	}

	bool operator < (const dr_star &other) {
		if (this->rating > other.rating) {
			return true;
		} else if (this->rating == other.rating) {
			if (this->name < other.name)
				return true;
		}
		return false;
	}

	dr_star & operator = (const float &other) {
		this->rating = other;
		this->name = "aaa\0";

		return *this;
	}
	friend std::ofstream &operator<<(std::ofstream &o_s, const dr_star &data);
};
std::ofstream &operator <<(std::ofstream &o_s, const dr_star &data) {
	o_s << data.name << ":";
	o_s << std::fixed << std::setprecision(3) <<data.rating;
	return o_s;
}

// type used for sorting by the number of rides the drivers fulfiled
struct dr_race {
	std::string name;
	int num_races;

	dr_race() {
		this->num_races = -1;
		this->name = "aaa\0";
	}

	dr_race(const int &num_races, const std::string &name) {
		this->name = name;
		this->num_races = num_races;
	}

	bool operator < (const dr_race &other) {
		if (this->num_races > other.num_races) {
			return true;
		} else if (this->num_races == other.num_races) {
			if (this->name < other.name)
				return true;
		}
		return false;
	}

	bool operator == (const dr_race &other) {
		if (this->name == other.name)
			return true;
		return false;
	}

	dr_race &operator += (const dr_race & other) {
		this->num_races = this->num_races + other.num_races;
		return *this;
	}

	dr_race & operator = (const int &other) {
		this->num_races = other;
		this->name = "aaa\0";

		return *this;
	}
	friend std::ofstream &operator<<(std::ofstream &o_r, const dr_race &data);
};

std::ofstream &operator <<(std::ofstream &o_r, const dr_race &data) {
	o_r << data.name << ":" << data.num_races;
	return o_r;
}

// type used for storing the information regarding the drivers
struct dr_status {
	std::string name, street;
	int dist = 0, num_races = 0;
	float num_stars = 0, rating = 0;
	bool is_online = 0;

	dr_status() { }
};

//////////////////////////////////////
// Skip_list

// the data found in one node of the skip_list
template <typename the_type> struct Skip_List_node {
	int idx = -1;
	int level = 0;
	the_type value;

    std::vector<Skip_List_node<the_type>*> level_point_next;
    std::vector<Skip_List_node<the_type>*> level_point_prev;
};

// defining the skip_list class
template <typename the_type> class Skip_List {
 private:
    // the borders of the skip_list
    Skip_List_node<the_type> *head;
    Skip_List_node<the_type> *tail;

    // elemnts regarding the structure of the skip_list
    int max_level;
    int elm_num;
    unsigned int seed = time(0);

 public:
    Skip_List() {
        max_level = 1;
        elm_num = 0;
        // constructor that creates and initializes the skip_list
        head = new Skip_List_node<the_type>;
        tail = new Skip_List_node<the_type>;

        head->level_point_next.push_back(tail);
        head->level_point_prev.push_back(nullptr);
        tail->level_point_next.push_back(nullptr);
        tail->level_point_prev.push_back(head);

        head->level = 1;
        tail->level = 1;

        head->value = -1;
        tail->value = -1;
    }
    // destructor of the skiplist
    ~Skip_List() {
    	Skip_List_node<the_type> *aux = head;
		Skip_List_node<the_type> *next_aux = head;
		do {
			next_aux = aux->level_point_next[0];
			delete aux;
			aux = next_aux;
		} while (next_aux != nullptr);
    }

    // function used for inserting elemnts in the skiplist
    void addinSkip_list(the_type data) {
    	++elm_num;
    	// updating the max possible level considering the number of elemts
    	if (elm_num  >= 1 << max_level) {
    		max_level++;

            head->level_point_next.push_back(tail);
            head->level_point_prev.push_back(nullptr);
            tail->level_point_next.push_back(nullptr);
            tail->level_point_prev.push_back(head);

            head->level = max_level;
            tail->level = max_level;
    	}
    	// creating a new node and fixing it's level
    	Skip_List_node<the_type> *aux = new Skip_List_node<the_type>;
    	Skip_List_node<the_type> *pos_finder = head;
    	int from_level = rand_r(&seed) % max_level;
    	int this_level = max_level;
    	aux->level_point_next.resize(from_level + 1);
    	aux->level_point_prev.resize(from_level + 1);
    	aux->level = from_level + 1;
    	aux->idx = 1;
    	aux->value = data;
    	// going through the skip_list and fining the right position where
    	// the new element should be
    	while (this_level > 0) {
    		--this_level;
    		bool ok = 1;
    		while (ok) {
    			if (pos_finder -> level_point_next[this_level]->idx == -1) {
    				ok = 0;
    			} else {
    				if (pos_finder -> level_point_next[this_level]->value <
    				data) {
    					pos_finder = pos_finder->level_point_next[this_level];
    				} else {
    					ok = 0;
    				}
    			}
    		}
    		// linking the new element with the other from the specific level
    		if (this_level <= from_level) {
    			aux->level_point_prev[this_level] = pos_finder;

    			aux->level_point_next[this_level] =
    			pos_finder -> level_point_next[this_level];

    			Skip_List_node<the_type> *temp;

    			temp = pos_finder -> level_point_next[this_level];
    			pos_finder -> level_point_next[this_level] = aux;

    			temp->level_point_prev[this_level] = aux;
    		}
    	}
    }

    // function used for finding a specific value in the skip_list
    Skip_List_node<the_type>* inSkip_list(the_type data) {
    	int this_level = max_level;
        Skip_List_node<the_type> *pos_finder = head;
        while (this_level > 0) {
            --this_level;
            bool ok = 1;
            while (ok) {
                if (pos_finder -> level_point_next[this_level]->idx == -1) {
                    ok = 0;
                } else {
                    if (pos_finder -> level_point_next[this_level]->value <
                    data) {
                        pos_finder = pos_finder->level_point_next[this_level];
                    } else {
                        ok = 0;
                    }
                }
            }
            if (pos_finder->level_point_next[this_level]->value == data)
            	return pos_finder->level_point_next[this_level];
        }
        // returning null_pointer if the value we were looking for is not to be
        // found
        return nullptr;
    }

    // funtion that removes elements from the skip_list
    void deletefromSkip_list(Skip_List_node<the_type>* pos_finder) {
    	Skip_List_node<the_type> *back_pos;
        Skip_List_node<the_type> *frnt_pos;
        --elm_num;

        // linking the remaining nodes
        for (int i = pos_finder->level - 1; i >= 0; --i) {
            back_pos = pos_finder->level_point_prev[i];
            frnt_pos = pos_finder->level_point_next[i];

            back_pos->level_point_next[i] = frnt_pos;
            frnt_pos->level_point_prev[i] = back_pos;
        }
        delete pos_finder;
        return;
    }

    // displaying the first elements of the skip_list, in the specific order
    void show_list(int how_many, std::ofstream& outf) {
    	Skip_List_node<the_type> *aux = head;
    	int i = 0;
    	// literally, how_many :)
    	while ( i < how_many ) {
    		aux = aux->level_point_next[0];
    		if (aux->level_point_next[0] == nullptr)
    			break;
    		outf << aux->value;
    		outf << " ";
    		++i;
    	}
    	outf << "\n";
    }

    // cheking if the value that is getting updated needs to be moved
    // used for keeping the elements in the correct order
    void chek_this_value(the_type data, the_type adding) {
    	Skip_List_node<the_type> *pos_finder = inSkip_list(data);
    	Skip_List_node<the_type> *adjecent;

    	pos_finder->value += adding;
    	the_type aux = pos_finder->value;

    	adjecent = pos_finder->level_point_prev[0];
    	if (adjecent->idx != -1) {
    		if (pos_finder->value < adjecent->value) {
    			deletefromSkip_list(pos_finder);
    			addinSkip_list(aux);
    			return;
    		}
    	}

    	adjecent = pos_finder->level_point_next[0];
    	if (adjecent->idx != -1) {
    		if (adjecent->value < pos_finder->value) {
    			deletefromSkip_list(pos_finder);
    			addinSkip_list(aux);
    			return;
    		}
    	}
    }
    // going to the next element
    Skip_List_node<the_type>* next(Skip_List_node<the_type>* node) {
    	return node->level_point_next[0];
    }
    // going to the previous element
    Skip_List_node<the_type>* prev(Skip_List_node<the_type>* node) {
    	return node->level_point_prev[0];
    }
    // setter for the size of the skip_list (number of elemts)
    int size() {
    	return elm_num;
    }
    // setters for the pointers to the front and the back of the skip_list
    Skip_List_node<the_type>* begin() {
    	return head->level_point_next[0];
    }
    Skip_List_node<the_type>* end() {
    	return tail;
    }
};

////////////////////////////////////
// Hash_Table

// the structure of the element fount in the hash_table
template <typename the_type> struct Node {
	the_type the_data;
	int node_number;

	Node() { }

	Node(the_type the_data, int node_number) {
		this->the_data = the_data;
		this->node_number = node_number;
	}
};

// absolute value
int abs_value(int x) {
	if (x < 0)
		return -x;
	return x;
}

// hashing function
template <typename the_type> unsigned int key_for_hash(the_type elemnet) {
	std::string aux_str(elemnet);
	unsigned int hash_value = 0;
	int p = 1;
	for (unsigned int i = 0; i < aux_str.length(); ++i) {
		hash_value += abs_value(aux_str[i] - 'Z') * p;
		p = p << 1;
	}
	hash_value %= prime_number;
	return hash_value;
}

// hashing class
template <typename the_data> class hashing_class {
	private:
		std::vector< std::list< Node< the_data > > > table_v;
		int elm_num = 0;

	public:
	~hashing_class() { }
	// adding an element to the hash_table and returning it's position
	int add_in_hash(the_data element) {
		unsigned int key = key_for_hash(element);
		if (table_v.size() <= key)
			table_v.resize(key + 1);
		table_v[key].push_back({element, elm_num});
		++elm_num;
		return elm_num - 1;
	}
	// looking for a specific element in the hash_table
	int return_poz(the_data element) {
		unsigned int key = key_for_hash(element);
		if (key < table_v.size()) {
			Node< the_data > p;
			typename std::list< Node< the_data > >::iterator it;
			it = table_v[key].begin();
			while (it != table_v[key].end()) {
				p = *it;
				if (p.the_data == element) {
					return p.node_number;
				}
				++it;
			}
		}
		// getting -1 if the value is not to be found
		return -1;
	}
	// setter for the number of elements in the hash_table
	int return_size() {
		return elm_num;
	}
};


/////////////////
// Solve class

class solver {
	private:
		int n, m, **distances;
		// the data I need to solve the tasks
		std::vector<dr_status> driver_list;
		std::queue<int> qu;
		hashing_class<std::string> H, dr_H;
		std::vector< Skip_List<int> > list_nodes;
		std::vector<std::string> street_names;

	public:
		// function used to read and process the initial map
		void read_data(std::ifstream& fin) {
			fin >> n >> m;
			int nod1, nod2;
			std::string str1, str2;
			distances = new int* [n];
			for (int i = 0; i < n; ++i) {
				distances[i] = new int[n];
			}
			// hashing all the streets in order to store them eficiently
			for (int i = 0; i < n; ++i) {
				fin >> str1;
				H.add_in_hash(str1);
				street_names.push_back(str1);
			}
			// linking the streets
			list_nodes.resize(n);
			for (int i = 0; i < m; ++i) {
				fin >> str1 >> str2;
				nod1 = H.return_poz(str1);
				nod2 = H.return_poz(str2);
				// using the skip_list in ordet to know the other streets that
				// are accessible from one street
				list_nodes[nod1].addinSkip_list(nod2);
			}
		}
		// breadth first search to find if another node can be accessed from
		// the starting position and, eventualy, the distance to it
		int B_F_S(const int &start, const int &finish) {
			int *dist = new int[n];
			std::fill(dist, dist + n, large_number);
			int nod, rez = 0;

			dist[start] = 0;
			qu.push(start);

			Skip_List_node<int>* impro_iter;

			while(!qu.empty()) {
				nod = qu.front();
				qu.pop();
				if (nod == finish)
					break;

				impro_iter = list_nodes[nod].begin();
				for (int i = 0; i < list_nodes[nod].size(); ++i) {
					if (dist[ impro_iter->value ] > 1 + dist[nod]) {
						dist[ impro_iter->value ] = 1 + dist[nod];
						qu.push(impro_iter->value);
					}
					impro_iter = list_nodes[nod].next(impro_iter);
				}
			}

			while(!qu.empty()) {
				qu.pop();
			}
			rez = dist[finish];
			delete [] dist;
			return rez;
		}
		// solver of task 1
		void task1_solver(std::ifstream& fin, std::ofstream& fout) {
			read_data(fin);
			int q1, nod2, nod1;
			std::string str1, str2;
			fin >> q1;
			for (int i = 0; i < q1; ++i) {
				fin >> str1 >> str2;
				nod1 = H.return_poz(str1);
				nod2 = H.return_poz(str2);
				if (large_number == B_F_S(nod1, nod2)) {
					fout << "n\n";
				} else {
					fout << "y\n";
				}
			}
		}
		// solver of thask 2 ... resembling to task 1
		void task2_solver(std::ifstream& fin, std::ofstream& fout) {
			int q2, nod2, nod1;
			std::string str1, str2;
			fin >> q2;
			for (int i = 0; i < q2; ++i) {
				fin >> str1 >> str2;
				nod1 = H.return_poz(str1);
				nod2 = H.return_poz(str2);
				int dist = B_F_S(nod1, nod2);
				if (large_number == dist) {
					fout << "-1\n";
				} else {
					fout << dist << "\n";
				}
			}
		}
		// function used to delete the linking between roads
		void cutting_road(Skip_List_node<int>* impro_iter, const int nod,
			Skip_List<int> &list_node) {
			Skip_List_node<int>* temp_del;
			while (impro_iter->value == nod) {
				temp_del = impro_iter;
				impro_iter = list_node.next(impro_iter);
				list_node.deletefromSkip_list(temp_del);
			}
		}
		// function used to compute all the possible distances between
		// ALL the nodes
		// mapping all the possible ways and all the possible ditances
		// and also keeping them in a proper order
		void total_distances() {
			Skip_List_node<int>* impro_iter;
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					distances[i][j] = large_number;
				}
				int nod = 0;
				distances[i][i] = 0;
				qu.push(i);
				// breadth first search from every single node to all others
				while (!qu.empty()) {
					nod = qu.front();
					qu.pop();
					impro_iter = list_nodes[nod].begin();

					for (int k = 0; k < list_nodes[nod].size(); ++k) {
						if (distances[i][ impro_iter->value ] >
							1 + distances[i][nod]) {
							distances[i][ impro_iter->value ] =
							1 + distances[i][nod];
							qu.push(impro_iter->value);
						}
						impro_iter = list_nodes[nod].next(impro_iter);
					}
				}
			}
		}
		// solver for task 3
		void task3_solver(std::ifstream& fin, std::ofstream& fout) {
			int q3, nod2, nod1, nod3, action, distsig, dist1, dist2;
			char c;
			std::string str1, str2, str3;
			Skip_List_node<int> *impro_iter, *impro_iter2;
			fin >> q3;
			for (int i = 0; i < q3; ++i) {
				fin >> c >> str1 >> str2 >> action;
				nod1 = H.return_poz(str1);
				nod2 = H.return_poz(str2);
				// fulfilling all the required actions
				if (c == 'c') {
					if (action == 0) {
						// making a link ; don't if there already is one
						impro_iter = list_nodes[nod1].inSkip_list(nod2);
						if (impro_iter == nullptr) {
							list_nodes[nod1].addinSkip_list(nod2);
						}
					} else if (action == 1) {
						// removing all the existent links between those two
						impro_iter = list_nodes[nod1].inSkip_list(nod2);
						if (impro_iter != nullptr) {
							cutting_road(impro_iter, nod2, list_nodes[nod1]);
						}
						impro_iter = list_nodes[nod2].inSkip_list(nod1);
						if (impro_iter != nullptr) {
							cutting_road(impro_iter, nod1, list_nodes[nod2]);
						}
					} else if (action == 2) {
						// road that goes both ways
						impro_iter = list_nodes[nod1].inSkip_list(nod2);
						if (impro_iter == nullptr) {
							list_nodes[nod1].addinSkip_list(nod2);
						}
						impro_iter = list_nodes[nod2].inSkip_list(nod1);
						if (impro_iter == nullptr) {
							list_nodes[nod2].addinSkip_list(nod1);
						}
					} else if (action == 3) {
						// flip the road's way :)
						impro_iter = list_nodes[nod1].inSkip_list(nod2);
						impro_iter2 = list_nodes[nod2].inSkip_list(nod1);
						if (impro_iter == nullptr && impro_iter2 != nullptr) {
							cutting_road(impro_iter2, nod1,
								list_nodes[nod2]);
							list_nodes[nod1].addinSkip_list(nod2);
						}
						if (impro_iter != nullptr && impro_iter2 == nullptr) {
							cutting_road(impro_iter, nod2,
								list_nodes[nod1]);
							list_nodes[nod2].addinSkip_list(nod1);
						}
					}
				} else {
					if (action == 0) {
						// checking if there is a way
						distsig =  B_F_S(nod1, nod2);
						if (large_number == distsig) {
							fout << "n\n";
						} else {
							fout << "y\n";
						}
					} else if (action == 1) {
						// checking the distance
						distsig = B_F_S(nod1, nod2);
						if (large_number == distsig) {
							fout << "-1\n";
						} else {
							fout << distsig << "\n";
						}
					} else if (action == 2) {
						// checking the distance passing by another location
						fin >> str3;
						nod3 = H.return_poz(str3);
						dist1 = B_F_S(nod1, nod3);
						dist2 = B_F_S(nod3, nod2);
						if (large_number == dist1 || large_number == dist2) {
							fout << "-1\n";
						} else {
							fout << dist1 + dist2 << "\n";
						}
					}
				}
			}
			// computing the distances
			// this will help reduce the running-time for the future tasks
			total_distances();
		}

		// introducing a new driver
		int make_new_driver(const std::string &the_name,
			std::vector<dr_status> &driver_list, Skip_List<dr_dist> &dist_list,
			Skip_List<dr_race> &race_list, Skip_List<dr_star> &star_list) {
			// inserting in the database
			int idx = dr_H.add_in_hash(the_name);
			driver_list.push_back({});
			driver_list[idx].name = the_name;

			// adding driver's data to the needed lists in order to keep them
			// sorted
			dist_list.addinSkip_list({0, the_name});
			race_list.addinSkip_list({0, the_name});
			star_list.addinSkip_list({0, 0, the_name});
			return idx;
		}

		// finding the suitable driver for this ride
		void get_short_race(int &mindist, const int &dr_this, const int &n1,
		int &dr_top, const std::vector<dr_status> &dr_list) {
			int n0 = H.return_poz(dr_list[dr_this].street);
			if (dr_top == -1) {
				dr_top = dr_this;
				mindist = distances[n0][n1];
				return;
			}
			if (mindist > distances[n0][n1]) {
				mindist = distances[n0][n1];
				dr_top = dr_this;
			} else {
				if (mindist == distances[n0][n1]) {
					if (dr_list[dr_this].rating > dr_list[dr_top].rating) {
						dr_top = dr_this;
					} else {
						if (dr_list[dr_this].rating ==
							dr_list[dr_top].rating) {
							if (dr_list[dr_this].name < dr_list[dr_top].name)
								dr_top = dr_this;
						}
					}
				}
			}
		}
		// the solution for task 4
		void task4_solver(std::ifstream& fin, std::ofstream& fout) {
			int q4, nod2, nod1, idx, mindist;
			Skip_List_node<int>* impro_iter;
			unsigned int num_cnt;
			std::string comm, the_name, str1, str2;
			float stars;
			// initializing the needed structures in order to solve this task
			Skip_List<dr_dist> dist_list;
			Skip_List<dr_race> race_list;
			Skip_List<dr_star> star_list;
			bool ok = 0;
			dr_dist dist_update;
			dr_race race_update;
			dr_star star_update;

			fin >> q4;

			for (int i = 0; i < q4; ++i) {
				fin >> comm;
				// getting all the commands
				if (comm == "d") {
					// driver gets online or joins in
					fin >> the_name >> str1;
					idx = dr_H.return_poz(the_name);
					if (idx == -1) {
						idx = make_new_driver(the_name, driver_list,
							dist_list, race_list, star_list);
					}
					driver_list[idx].is_online = 1;
					driver_list[idx].street = str1;
				} else if (comm == "b") {
					// driver gets offline
					fin >> the_name;
					idx = dr_H.return_poz(the_name);
					if (idx == -1) {
						idx = make_new_driver(the_name, driver_list,
							dist_list, race_list, star_list);
					}
					driver_list[idx].is_online = 0;
				} else if (comm == "r") {
					// notification for a new ride
					fin >> str1 >> str2 >> stars;
					nod1 = H.return_poz(str1);
					nod2 = H.return_poz(str2);
					mindist = large_number + 10;
					idx = -1;

					// looking for the suitable driver
					for (unsigned int j = 0; j < driver_list.size() ; ++j) {
						if (driver_list[j].is_online == 1) {
							get_short_race(mindist, j, nod1, idx,
							driver_list);
						}
					}
					// text when no driver are available to pick the ride
					if (mindist >= large_number) {
						fout << "Soferi indisponibili\n";
						continue;
					}
					ok = 1;
					if (distances[nod1][nod2] == large_number) {
						impro_iter = list_nodes[nod2].begin();
						ok = 0;
						// if the finishing point is not accessible
						// try finding a close spot
						for (int k = 0; k < list_nodes[nod2].size();
							++k) {
							if (distances[nod1][ impro_iter->value ] <
								large_number) {
								nod2 = impro_iter->value;
								str2 = street_names[nod2];
								ok = 1;
								break;
							} else {
								impro_iter = list_nodes[nod2].next(impro_iter);
							}
						}
					}
					// text when the destination woun't be reached
					if (ok == 0) {
						fout << "Destinatie inaccesibila\n";
						continue;
					}

					// updating all the lists
					dist_list.chek_this_value({driver_list[idx].dist,
						driver_list[idx].name},
						{mindist + distances[nod1][nod2],
							driver_list[idx].name});

					race_list.chek_this_value({driver_list[idx].num_races,
						driver_list[idx].name}, {1, driver_list[idx].name});

					star_list.chek_this_value({driver_list[idx].num_races,
						driver_list[idx].num_stars, driver_list[idx].name},
						{1, stars, driver_list[idx].name});

					// updating the data base
					driver_list[idx].street = str2;
					driver_list[idx].dist += mindist + distances[nod1][nod2];
					++driver_list[idx].num_races;
					driver_list[idx].num_stars += stars;
					driver_list[idx].rating =
					driver_list[idx].num_stars /
					driver_list[idx].num_races;

					// printing the data
				} else if (comm == "top_rating") {
					fin >> num_cnt;
					star_list.show_list(num_cnt, fout);
				} else if (comm == "top_dist") {
					fin >> num_cnt;
					dist_list.show_list(num_cnt, fout);
				} else if (comm == "top_rides") {
					fin >> num_cnt;
					race_list.show_list(num_cnt, fout);
				} else if (comm == "info") {
					fin >> the_name;
					idx = dr_H.return_poz(the_name);

					fout << the_name << ": " << driver_list[idx].street << " "
					<< driver_list[idx].rating << " "
					<< driver_list[idx].num_races << " "
					<< driver_list[idx].dist;
					if (driver_list[idx].is_online == 1)
						fout << " online\n";
					else
						fout << " ofline\n";
				}
			}
		}

		// solver for the task 5
		void task5_solver(std::ifstream& fin, std::ofstream& fout) {
			int fuel, idx, no_str, node;
			std::string name;
			bool *ok_way = new bool[n]();

			// processing the input data
			fin >> fuel >> name >> no_str;
			idx = dr_H.return_poz(name);
			node = H.return_poz(driver_list[idx].street);

			// marking the spots of interest for the driver
			for (int i = 0; i < no_str; ++i) {
				fin >> name;
				ok_way[H.return_poz(name)] = 1;
			}

			if (ok_way[node] == 1) {
				fout << street_names[node] << " ";
			}

			// checking the map of the distances to see the ceratin spots
			// that can be reached
			for (int i = 1; i <= fuel; ++i) {
				for (int j = 0; j < n; ++j) {
					if (distances[node][j] == i && ok_way[j] == 1)
						fout << street_names[j] << " ";
				}
			}

			for (int i = 0; i < n; ++i) {
				delete [] distances[i];
			}
			delete [] ok_way;
			delete [] distances;
		}
};

#endif  // SOLVER_H_
