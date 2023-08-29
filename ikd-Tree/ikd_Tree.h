#pragma once
#include <stdio.h>
#include <queue>
#include <pthread.h>
#include <chrono>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <algorithm>
#include <memory>
#include <pcl/point_types.h>

#define EPSS 1e-6
#define Minimal_Unbalanced_Tree_Size 10
#define Multi_Thread_Rebuild_Point_Num 1500
#define DOWNSAMPLE_SWITCH true
#define ForceRebuildPercentage 0.2
#define Q_LEN 1000000

using namespace std;

struct ikdTree_PointType
{
    float x, y, z;
    ikdTree_PointType(float px = 0.0f, float py = 0.0f, float pz = 0.0f)
    {
        x = px;
        y = py;
        z = pz;
    }
};

struct BoxPointType
{
    float vertex_min[3];
    float vertex_max[3];
};

enum operation_set
{
    ADD_POINT,
    DELETE_POINT,
    DELETE_BOX,
    ADD_BOX,
    DOWNSAMPLE_DELETE,
    PUSH_DOWN
};

enum delete_point_storage_set
{
    NOT_RECORD,
    DELETE_POINTS_REC,
    MULTI_THREAD_REC
};

template <typename T>
class MANUAL_Q
{
private:
    int head = 0, tail = 0, counter = 0;
    T q[Q_LEN];
    bool is_empty;

public:
    void pop();
    T front();
    T back();
    void clear();
    void push(T op);
    bool empty();
    int size();
};

template <typename PointType>
class KD_TREE
{
public:
    using PointVector = vector<PointType, Eigen::aligned_allocator<PointType>>;
    using Ptr = shared_ptr<KD_TREE<PointType>>;
    struct KD_TREE_NODE
    {
        PointType point;           // 数据点
        uint8_t division_axis;     // 分割轴
        int TreeSize = 1;          // 总节点数
        int invalid_point_num = 0; // label 为删除的点的数量
        int down_del_num = 0;
        bool point_deleted = false; // 当前节点是否标记为删除
        bool tree_deleted = false;  // 整个(子)树是否标记为删除
        bool point_downsample_deleted = false;
        bool tree_downsample_deleted = false;
        bool need_push_down_to_left = false;
        bool need_push_down_to_right = false;
        bool working_flag = false;
        float radius_sq;
        pthread_mutex_t push_down_mutex_lock;
        float node_range_x[2], node_range_y[2], node_range_z[2]; // tree 对应的包络Boxs
        KD_TREE_NODE *left_son_ptr = nullptr;                    // 左子树
        KD_TREE_NODE *right_son_ptr = nullptr;                   // 右子树
        KD_TREE_NODE *father_ptr = nullptr;                      // 父子树
        // For paper data record
        float alpha_del;
        float alpha_bal;
    };

    struct Operation_Logger_Type
    {
        PointType point;
        BoxPointType boxpoint;
        bool tree_deleted, tree_downsample_deleted;
        operation_set op;
    };

    /**
     * @brief 堆节点的类型
     *
     */
    struct PointType_CMP
    {
        PointType point;
        float dist = 0.0;
        PointType_CMP(PointType p = PointType(), float d = INFINITY)
        {
            this->point = p;
            this->dist = d;
        };
        bool operator<(const PointType_CMP &a) const
        {
            if (fabs(dist - a.dist) < 1e-10)
                return point.x < a.point.x;
            else
                return dist < a.dist;
        }
    };

    /**
     * @brief 大根堆
     *
     */
    class MANUAL_HEAP
    {
    public:
        /**
         * @brief 初始化一个大根堆
         *
         * @param max_capacity  设定堆的最大容量
         */
        MANUAL_HEAP(int max_capacity = 100)
        {
            cap = max_capacity;
            heap = new PointType_CMP[max_capacity];
            heap_size = 0;
        }

        /**
         * @brief 释放兑数组的内存
         *
         */
        ~MANUAL_HEAP() { delete[] heap; }

        /**
         * @brief 从大根堆中弹出（删除）堆顶元素
         *
         */
        void pop()
        {
            if (heap_size == 0)
                return;
            heap[0] = heap[heap_size - 1];
            heap_size--;
            MoveDown(0);
            return;
        }

        /**
         * @brief 获取大根堆的最大值
         *
         * @return PointType_CMP
         */
        PointType_CMP top() { return heap[0]; }

        /**
         * @brief 向大根堆中插入新元素，如果堆已满，无法添加新元素
         *
         * @param point
         */
        void push(PointType_CMP point)
        {
            if (heap_size >= cap)
                return;
            heap[heap_size] = point;
            FloatUp(heap_size);
            heap_size++;
            return;
        }

        /**
         * @brief 获取当前堆中的节点数目
         *
         * @return int
         */
        int size() { return heap_size; }

        /**
         * @brief 清除整个堆
         *
         */
        void clear() { heap_size = 0; }

    private:
        int heap_size = 0;   // 堆中当前元素的数量
        int cap = 0;         // 堆的最大容量
        PointType_CMP *heap; // 用于存储堆元素的数组

        /**
         * @brief 大根堆中的节点下沉操作
         *
         * @param heap_index  要下沉的堆节点index
         */
        void MoveDown(int heap_index)
        {
            int l = heap_index * 2 + 1;           // 计算左子节点的索引
            PointType_CMP tmp = heap[heap_index]; // 保存当前节点的值

            // 不断调整堆的结构
            while (l < heap_size)
            {
                // 检查是否存在右子节点，并且右子节点的值大于左子节点，如果右子节点更大，则选择右子节点
                if (l + 1 < heap_size && heap[l] < heap[l + 1])
                    l++;

                // 比较当前节点的值和选定的子节点的值，如果当前节点的值小于子节点的值
                if (tmp < heap[l])
                {
                    heap[heap_index] = heap[l]; // 将子节点上移，覆盖当前节点
                    heap_index = l;             // 更新当前节点索引
                    l = heap_index * 2 + 1;     // 计算左子节点索引
                }
                else
                    break;
            }
            heap[heap_index] = tmp;
            return;
        }

        /**
         * @brief 大根堆中的节点上浮操作
         *
         * @param heap_index 要上浮的堆节点index
         */
        void FloatUp(int heap_index)
        {
            int ancestor = (heap_index - 1) / 2;  // 计算父节点的索引
            PointType_CMP tmp = heap[heap_index]; // 保存当前节点的值

            // 不断调整大根堆的结构
            while (heap_index > 0)
            {
                // 比较当前节点的值和父节点的值，如果父节点的值小于当前节点的值
                if (heap[ancestor] < tmp)
                {
                    heap[heap_index] = heap[ancestor]; // 将父节点下移，覆盖当前节点
                    heap_index = ancestor;             // 更新当前节点索引
                    ancestor = (heap_index - 1) / 2;   // 计算新的父节点索引
                }
                else
                    break;
            }
            heap[heap_index] = tmp; // 将当前节点放置在最终位置上
            return;
        }
    };

private:
    // Multi-thread Tree Rebuild
    bool termination_flag = false;
    bool rebuild_flag = false;
    pthread_t rebuild_thread;
    pthread_mutex_t termination_flag_mutex_lock, rebuild_ptr_mutex_lock, working_flag_mutex, search_flag_mutex;
    pthread_mutex_t rebuild_logger_mutex_lock, points_deleted_rebuild_mutex_lock;
    // queue<Operation_Logger_Type> Rebuild_Logger;
    MANUAL_Q<Operation_Logger_Type> Rebuild_Logger;
    PointVector Rebuild_PCL_Storage;
    KD_TREE_NODE **Rebuild_Ptr = nullptr;
    int search_mutex_counter = 0;
    static void *multi_thread_ptr(void *arg);
    void multi_thread_rebuild();
    void start_thread();
    void stop_thread();
    void run_operation(KD_TREE_NODE **root, Operation_Logger_Type operation);
    // KD Tree Functions and augmented variables
    int Treesize_tmp = 0, Validnum_tmp = 0;
    float alpha_bal_tmp = 0.5, alpha_del_tmp = 0.0;
    float delete_criterion_param = 0.5f;  // 删除剪枝
    float balance_criterion_param = 0.7f; // 平衡剪枝
    float downsample_size = 0.2f;         // 下采样分辨率
    bool Delete_Storage_Disabled = false;
    KD_TREE_NODE *STATIC_ROOT_NODE = nullptr;
    PointVector Points_deleted;
    PointVector Downsample_Storage;
    PointVector Multithread_Points_deleted;
    void InitTreeNode(KD_TREE_NODE *root);
    void Test_Lock_States(KD_TREE_NODE *root);
    void BuildTree(KD_TREE_NODE **root, int l, int r, PointVector &Storage);
    void Rebuild(KD_TREE_NODE **root);
    int Delete_by_range(KD_TREE_NODE **root, BoxPointType boxpoint, bool allow_rebuild, bool is_downsample);
    void Delete_by_point(KD_TREE_NODE **root, PointType point, bool allow_rebuild);
    void Add_by_point(KD_TREE_NODE **root, PointType point, bool allow_rebuild, int father_axis);
    void Add_by_range(KD_TREE_NODE **root, BoxPointType boxpoint, bool allow_rebuild);
    void Search(KD_TREE_NODE *root, int k_nearest, PointType point, MANUAL_HEAP &q, double max_dist); // priority_queue<PointType_CMP>
    void Search_by_range(KD_TREE_NODE *root, BoxPointType boxpoint, PointVector &Storage);
    void Search_by_radius(KD_TREE_NODE *root, PointType point, float radius, PointVector &Storage);
    bool Criterion_Check(KD_TREE_NODE *root);
    void Push_Down(KD_TREE_NODE *root);
    void Update(KD_TREE_NODE *root);
    void delete_tree_nodes(KD_TREE_NODE **root);
    void downsample(KD_TREE_NODE **root);
    bool same_point(PointType a, PointType b);
    float calc_dist(PointType a, PointType b);
    float calc_box_dist(KD_TREE_NODE *node, PointType point);
    static bool point_cmp_x(PointType a, PointType b);
    static bool point_cmp_y(PointType a, PointType b);
    static bool point_cmp_z(PointType a, PointType b);

public:
    KD_TREE(float delete_param = 0.5, float balance_param = 0.6, float box_length = 0.2);
    ~KD_TREE();
    void Set_delete_criterion_param(float delete_param);
    void Set_balance_criterion_param(float balance_param);
    void set_downsample_param(float box_length);
    void InitializeKDTree(float delete_param = 0.5, float balance_param = 0.7, float box_length = 0.2);
    int size();
    int validnum();
    void root_alpha(float &alpha_bal, float &alpha_del);
    void Build(PointVector point_cloud);
    void Nearest_Search(PointType point, int k_nearest, PointVector &Nearest_Points, vector<float> &Point_Distance, double max_dist = INFINITY);
    void Box_Search(const BoxPointType &Box_of_Point, PointVector &Storage);
    void Radius_Search(PointType point, const float radius, PointVector &Storage);
    int Add_Points(PointVector &PointToAdd, bool downsample_on);
    void Add_Point_Boxes(vector<BoxPointType> &BoxPoints);
    void Delete_Points(PointVector &PointToDel);
    int Delete_Point_Boxes(vector<BoxPointType> &BoxPoints);
    void flatten(KD_TREE_NODE *root, PointVector &Storage, delete_point_storage_set storage_type);
    void acquire_removed_points(PointVector &removed_points);
    BoxPointType tree_range();
    PointVector PCL_Storage;
    KD_TREE_NODE *Root_Node = nullptr;
    int max_queue_size = 0;
};
