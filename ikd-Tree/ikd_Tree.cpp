#include "ikd_Tree.h"

/*
Description: ikd-Tree: an incremental k-d tree for robotic applications
Author: Yixi Cai
email: yixicai@connect.hku.hk
*/

/**
 * @brief Construct a new kd tree<pointtype>::kd tree object
 *
 * @tparam PointType
 * @param delete_param  删除标准
 * @param balance_param 平衡标准
 * @param box_length    下采样分辨率
 */
template <typename PointType>
KD_TREE<PointType>::KD_TREE(float delete_param, float balance_param, float box_length)
{
    delete_criterion_param = delete_param;
    balance_criterion_param = balance_param;
    downsample_size = box_length;
    Rebuild_Logger.clear();
    termination_flag = false;
    start_thread();
}

/**
 * @brief Destroy the kd tree<pointtype>::kd tree object
 *
 * @tparam PointType
 */
template <typename PointType>
KD_TREE<PointType>::~KD_TREE()
{
    stop_thread();
    Delete_Storage_Disabled = true;
    delete_tree_nodes(&Root_Node);
    PointVector().swap(PCL_Storage);
    Rebuild_Logger.clear();
}

/**
 * @brief 设置 delete_criterion_param
 *
 * @tparam PointType
 * @param delete_param  删除标准
 */
template <typename PointType>
void KD_TREE<PointType>::Set_delete_criterion_param(float delete_param)
{
    delete_criterion_param = delete_param;
}

/**
 * @brief 设置 balance_criterion_param
 *
 * @tparam PointType
 * @param balance_param 平衡标准
 */
template <typename PointType>
void KD_TREE<PointType>::Set_balance_criterion_param(float balance_param)
{
    balance_criterion_param = balance_param;
}

/**
 * @brief 设置下采样分辨率
 *
 * @tparam PointType
 * @param downsample_param 下采样分辨率
 */
template <typename PointType>
void KD_TREE<PointType>::set_downsample_param(float downsample_param)
{
    downsample_size = downsample_param;
}

/**
 * @brief 初始化KDTree的一些基本参数
 *
 * @tparam PointType
 * @param delete_param  删除标准
 * @param balance_param 平衡标准
 * @param box_length    下采样分辨率
 */
template <typename PointType>
void KD_TREE<PointType>::InitializeKDTree(float delete_param, float balance_param, float box_length)
{
    Set_delete_criterion_param(delete_param);
    Set_balance_criterion_param(balance_param);
    set_downsample_param(box_length);
}

/**
 * @brief 初始化一个KDTree的节点
 *
 * @tparam PointType
 * @param root KDTree节点指针
 */
template <typename PointType>
void KD_TREE<PointType>::InitTreeNode(KD_TREE_NODE *root)
{
    root->point.x = 0.0f;
    root->point.y = 0.0f;
    root->point.z = 0.0f;
    root->node_range_x[0] = 0.0f;
    root->node_range_x[1] = 0.0f;
    root->node_range_y[0] = 0.0f;
    root->node_range_y[1] = 0.0f;
    root->node_range_z[0] = 0.0f;
    root->node_range_z[1] = 0.0f;
    root->division_axis = 0;
    root->father_ptr = nullptr;
    root->left_son_ptr = nullptr;
    root->right_son_ptr = nullptr;
    root->TreeSize = 0;
    root->invalid_point_num = 0;
    root->down_del_num = 0;
    root->point_deleted = false;
    root->tree_deleted = false;
    root->need_push_down_to_left = false;
    root->need_push_down_to_right = false;
    root->point_downsample_deleted = false;
    root->working_flag = false;
    pthread_mutex_init(&(root->push_down_mutex_lock), NULL);
}

/**
 * @brief   获取 KD 树的节点数量。
 *          它可能直接从 Root_Node 获取节点数量，也可能在重建操作期间获取 Treesize_tmp 的值
 *
 * @tparam PointType
 * @return int
 */
template <typename PointType>
int KD_TREE<PointType>::size()
{
    int s = 0;
    if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node)
    {
        if (Root_Node != nullptr)
        {
            return Root_Node->TreeSize;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        if (!pthread_mutex_trylock(&working_flag_mutex))
        {
            s = Root_Node->TreeSize;
            pthread_mutex_unlock(&working_flag_mutex);
            return s;
        }
        else
        {
            return Treesize_tmp;
        }
    }
}

/**
 * @brief 获取整个 KD 树的包络范围
 *
 * @tparam PointType
 * @return BoxPointType
 */
template <typename PointType>
BoxPointType KD_TREE<PointType>::tree_range()
{
    BoxPointType range;

    // 检查当前是否正在进行重建操作，以及重建操作是否涉及根节点
    if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node)
    {
        // 如果根节点存在，将 range 的边界设置为根节点的包络范围。
        if (Root_Node != nullptr)
        {
            range.vertex_min[0] = Root_Node->node_range_x[0];
            range.vertex_min[1] = Root_Node->node_range_y[0];
            range.vertex_min[2] = Root_Node->node_range_z[0];
            range.vertex_max[0] = Root_Node->node_range_x[1];
            range.vertex_max[1] = Root_Node->node_range_y[1];
            range.vertex_max[2] = Root_Node->node_range_z[1];
        }
        // 如果根节点为空，使用 memset 函数将 range 的内存全部置为零。
        else
        {
            memset(&range, 0, sizeof(range));
        }
    }
    // 如果 Rebuild_Ptr 不为空且指向根节点
    else
    {
        // 使用互斥锁 working_flag_mutex 来保护多线程访问。
        if (!pthread_mutex_trylock(&working_flag_mutex))
        {
            // 将 range 的边界设置为根节点的包络范围。
            range.vertex_min[0] = Root_Node->node_range_x[0];
            range.vertex_min[1] = Root_Node->node_range_y[0];
            range.vertex_min[2] = Root_Node->node_range_z[0];
            range.vertex_max[0] = Root_Node->node_range_x[1];
            range.vertex_max[1] = Root_Node->node_range_y[1];
            range.vertex_max[2] = Root_Node->node_range_z[1];
            pthread_mutex_unlock(&working_flag_mutex);
        }
        else
        {
            memset(&range, 0, sizeof(range));
        }
    }
    return range;
}

/**
 * @brief 获取有效节点的数量，即不被标记为删除的节点数目
 *
 * @tparam PointType
 * @return int
 */
template <typename PointType>
int KD_TREE<PointType>::validnum()
{
    int s = 0;
    // 检查当前是否正在进行重建操作，以及重建操作是否涉及根节点
    if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node)
    {
        if (Root_Node != nullptr)
            return (Root_Node->TreeSize - Root_Node->invalid_point_num);
        else
            return 0;
    }
    else
    {
        if (!pthread_mutex_trylock(&working_flag_mutex))
        {
            s = Root_Node->TreeSize - Root_Node->invalid_point_num;
            pthread_mutex_unlock(&working_flag_mutex);
            return s;
        }
        else
        {
            return -1;
        }
    }
}

/**
 * @brief 获取平衡重建触发因子和删除重建触发因子
 *
 * @tparam PointType
 * @param alpha_bal 平衡重建触发因子
 * @param alpha_del 删除重建触发因子
 */
template <typename PointType>
void KD_TREE<PointType>::root_alpha(float &alpha_bal, float &alpha_del)
{
    // 没有正在进行第二线程重建
    if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node)
    {
        alpha_bal = Root_Node->alpha_bal;
        alpha_del = Root_Node->alpha_del;
        return;
    }
    // 正在进行第二线程重建
    else
    {
        // 上锁成功
        if (!pthread_mutex_trylock(&working_flag_mutex))
        {
            alpha_bal = Root_Node->alpha_bal;
            alpha_del = Root_Node->alpha_del;
            pthread_mutex_unlock(&working_flag_mutex);
            return;
        }
        // 上锁失败，那重建之前备份的旧的两个因子数值
        else
        {
            alpha_bal = alpha_bal_tmp;
            alpha_del = alpha_del_tmp;
            return;
        }
    }
}

/**
 * @brief 用于初始化互斥锁并创建一个新线程，以执行多线程重建操作
 *
 * @tparam PointType
 */
template <typename PointType>
void KD_TREE<PointType>::start_thread()
{
    pthread_mutex_init(&termination_flag_mutex_lock, NULL);
    pthread_mutex_init(&rebuild_ptr_mutex_lock, NULL);
    pthread_mutex_init(&rebuild_logger_mutex_lock, NULL);
    pthread_mutex_init(&points_deleted_rebuild_mutex_lock, NULL);
    pthread_mutex_init(&working_flag_mutex, NULL);
    pthread_mutex_init(&search_flag_mutex, NULL);
    pthread_create(&rebuild_thread, NULL, multi_thread_ptr, (void *)this);
    printf("Multi thread started \n");
}

/**
 * @brief 停止多线程重建操作。它会设置终止标志，等待重建线程结束，然后销毁相关的互斥锁。这确保了在停止多线程操作时的安全资源释放。
 *
 * @tparam PointType
 */
template <typename PointType>
void KD_TREE<PointType>::stop_thread()
{
    pthread_mutex_lock(&termination_flag_mutex_lock);
    termination_flag = true;
    pthread_mutex_unlock(&termination_flag_mutex_lock);
    if (rebuild_thread)
        pthread_join(rebuild_thread, NULL);
    pthread_mutex_destroy(&termination_flag_mutex_lock);
    pthread_mutex_destroy(&rebuild_logger_mutex_lock);
    pthread_mutex_destroy(&rebuild_ptr_mutex_lock);
    pthread_mutex_destroy(&points_deleted_rebuild_mutex_lock);
    pthread_mutex_destroy(&working_flag_mutex);
    pthread_mutex_destroy(&search_flag_mutex);
}

/**
 * @brief 另一个线程要执行的函数
 *
 * @tparam PointType
 * @param arg 传入参数
 * @return void*
 */
template <typename PointType>
void *KD_TREE<PointType>::multi_thread_ptr(void *arg)
{
    KD_TREE *handle = (KD_TREE *)arg;
    handle->multi_thread_rebuild();
    return nullptr;
}

/**
 * @brief 实现了一个多线程的 KD 树重建过程，保证了重建过程中各个操作的同步和正确性。
 *
 * @tparam PointType
 */
template <typename PointType>
void KD_TREE<PointType>::multi_thread_rebuild()
{
    bool terminated = false;
    KD_TREE_NODE *father_ptr, **new_node_ptr;
    pthread_mutex_lock(&termination_flag_mutex_lock);
    terminated = termination_flag; // 判断是否已经停止多线程
    pthread_mutex_unlock(&termination_flag_mutex_lock);

    // 首先通过 terminated 来判断是否终止多线程。如果没有终止，则进入多线程重建流程。
    while (!terminated)
    {
        pthread_mutex_lock(&rebuild_ptr_mutex_lock);
        pthread_mutex_lock(&working_flag_mutex);

        // Rebuild_Ptr不为空，需要第二线程重建
        if (Rebuild_Ptr != nullptr)
        {
            /* Traverse and copy */
            if (!Rebuild_Logger.empty())
            {
                printf("\n\n\n\n\n\n\n\n\n\n\n ERROR!!! \n\n\n\n\n\n\n\n\n");
            }

            // 将 rebuild_flag 设置为 true，以指示重建正在进行
            rebuild_flag = true;

            // 备份当前重建指针指向的节点信息(重建前的信息)，包括树的大小、有效点数等。
            if (*Rebuild_Ptr == Root_Node)
            {
                Treesize_tmp = Root_Node->TreeSize;
                Validnum_tmp = Root_Node->TreeSize - Root_Node->invalid_point_num;
                alpha_bal_tmp = Root_Node->alpha_bal;
                alpha_del_tmp = Root_Node->alpha_del;
            }
            KD_TREE_NODE *old_root_node = (*Rebuild_Ptr);
            father_ptr = (*Rebuild_Ptr)->father_ptr;
            PointVector().swap(Rebuild_PCL_Storage);

            // Lock Search
            // search_mutex_counter != 0 ,主线程正在使用搜索，继续等待，否则将 search_mutex_counter 置为 -1, 表示多线程重建占用
            pthread_mutex_lock(&search_flag_mutex);
            while (search_mutex_counter != 0) // 主线程正在使用搜索
            {
                pthread_mutex_unlock(&search_flag_mutex);
                usleep(1);
                pthread_mutex_lock(&search_flag_mutex);
            }
            search_mutex_counter = -1;
            pthread_mutex_unlock(&search_flag_mutex);

            // Lock deleted points cache 锁住删除节点统计变量 Multithread_Points_deleted 和 Points_deleted
            pthread_mutex_lock(&points_deleted_rebuild_mutex_lock);
            flatten(*Rebuild_Ptr, Rebuild_PCL_Storage, MULTI_THREAD_REC);
            // Unlock deleted points cache
            pthread_mutex_unlock(&points_deleted_rebuild_mutex_lock);

            // Unlock Search
            // 将 search_mutex_counter 置为 0, 表示多线程重建不再用
            pthread_mutex_lock(&search_flag_mutex);
            search_mutex_counter = 0;
            pthread_mutex_unlock(&search_flag_mutex);
            pthread_mutex_unlock(&working_flag_mutex);

            /* Rebuild and update missed operations*/
            Operation_Logger_Type Operation;
            KD_TREE_NODE *new_root_node = nullptr;
            if (int(Rebuild_PCL_Storage.size()) > 0)
            {
                // 先重建树
                BuildTree(&new_root_node, 0, Rebuild_PCL_Storage.size() - 1, Rebuild_PCL_Storage);

                // 后将重建logger中的操作补上
                // Rebuild has been done. Updates the blocked operations into the new tree
                pthread_mutex_lock(&working_flag_mutex);
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                int tmp_counter = 0;
                while (!Rebuild_Logger.empty())
                {
                    Operation = Rebuild_Logger.front();
                    max_queue_size = max(max_queue_size, Rebuild_Logger.size());
                    Rebuild_Logger.pop();
                    pthread_mutex_unlock(&rebuild_logger_mutex_lock);
                    pthread_mutex_unlock(&working_flag_mutex);
                    run_operation(&new_root_node, Operation);
                    tmp_counter++;
                    if (tmp_counter % 10 == 0)
                        usleep(1);
                    pthread_mutex_lock(&working_flag_mutex);
                    pthread_mutex_lock(&rebuild_logger_mutex_lock);
                }
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }

            // 将新构建的树替换原来的树
            /* Replace to original tree*/
            // pthread_mutex_lock(&working_flag_mutex);
            pthread_mutex_lock(&search_flag_mutex);
            while (search_mutex_counter != 0) // 主线程正在使用搜索
            {
                pthread_mutex_unlock(&search_flag_mutex);
                usleep(1);
                pthread_mutex_lock(&search_flag_mutex);
            }
            search_mutex_counter = -1;
            pthread_mutex_unlock(&search_flag_mutex);
            if (father_ptr->left_son_ptr == *Rebuild_Ptr)
            {
                father_ptr->left_son_ptr = new_root_node;
            }
            else if (father_ptr->right_son_ptr == *Rebuild_Ptr)
            {
                father_ptr->right_son_ptr = new_root_node;
            }
            else
            {
                throw "Error: Father ptr incompatible with current node\n";
            }
            if (new_root_node != nullptr)
                new_root_node->father_ptr = father_ptr;
            (*Rebuild_Ptr) = new_root_node;

            int valid_old = old_root_node->TreeSize - old_root_node->invalid_point_num;
            int valid_new = new_root_node->TreeSize - new_root_node->invalid_point_num;
            if (father_ptr == STATIC_ROOT_NODE)
                Root_Node = STATIC_ROOT_NODE->left_son_ptr;
            KD_TREE_NODE *update_root = *Rebuild_Ptr;
            while (update_root != nullptr && update_root != Root_Node)
            {
                update_root = update_root->father_ptr;
                if (update_root->working_flag)
                    break;
                if (update_root == update_root->father_ptr->left_son_ptr && update_root->father_ptr->need_push_down_to_left)
                    break;
                if (update_root == update_root->father_ptr->right_son_ptr && update_root->father_ptr->need_push_down_to_right)
                    break;
                Update(update_root);
            }
            pthread_mutex_lock(&search_flag_mutex);
            search_mutex_counter = 0;
            pthread_mutex_unlock(&search_flag_mutex);
            Rebuild_Ptr = nullptr;
            pthread_mutex_unlock(&working_flag_mutex);
            rebuild_flag = false;
            /* Delete discarded tree nodes */
            delete_tree_nodes(&old_root_node);
        }
        // 不需要重建
        else
        {
            pthread_mutex_unlock(&working_flag_mutex);
        }
        pthread_mutex_unlock(&rebuild_ptr_mutex_lock);
        pthread_mutex_lock(&termination_flag_mutex_lock);
        terminated = termination_flag;
        pthread_mutex_unlock(&termination_flag_mutex_lock);
        usleep(100);
    }
    printf("Rebuild thread terminated normally\n");
}

/**
 * @brief 根据传入的操作类型，选择性地调用不同的函数来执行相应的操作，这些操作涉及树的结构修改、添加、删除等
 *
 * @tparam PointType
 * @param root  当前(子)树根节点
 * @param operation 操作类型
 */
template <typename PointType>
void KD_TREE<PointType>::run_operation(KD_TREE_NODE **root, Operation_Logger_Type operation)
{
    // 所有的操作都不允许重建
    switch (operation.op)
    {
    case ADD_POINT:
        Add_by_point(root, operation.point, false, (*root)->division_axis);
        break;
    case ADD_BOX:
        Add_by_range(root, operation.boxpoint, false);
        break;
    case DELETE_POINT:
        Delete_by_point(root, operation.point, false);
        break;
    case DELETE_BOX:
        Delete_by_range(root, operation.boxpoint, false, false);
        break;
    case DOWNSAMPLE_DELETE:
        Delete_by_range(root, operation.boxpoint, false, true);
        break;
    case PUSH_DOWN:
        (*root)->tree_downsample_deleted |= operation.tree_downsample_deleted;
        (*root)->point_downsample_deleted |= operation.tree_downsample_deleted;
        (*root)->tree_deleted = operation.tree_deleted || (*root)->tree_downsample_deleted;
        (*root)->point_deleted = (*root)->tree_deleted || (*root)->point_downsample_deleted;
        if (operation.tree_downsample_deleted)
            (*root)->down_del_num = (*root)->TreeSize;
        if (operation.tree_deleted)
            (*root)->invalid_point_num = (*root)->TreeSize;
        else
            (*root)->invalid_point_num = (*root)->down_del_num;
        (*root)->need_push_down_to_left = true;
        (*root)->need_push_down_to_right = true;
        break;
    default:
        break;
    }
}

/**
 * @brief 根据给定点集合构建一棵KDTree
 *
 * @tparam PointType  给定点集合
 * @param point_cloud
 */
template <typename PointType>
void KD_TREE<PointType>::Build(PointVector point_cloud)
{
    if (Root_Node != nullptr)
    {
        delete_tree_nodes(&Root_Node);
    }
    if (point_cloud.size() == 0)
        return;
    STATIC_ROOT_NODE = new KD_TREE_NODE;
    InitTreeNode(STATIC_ROOT_NODE);
    BuildTree(&STATIC_ROOT_NODE->left_son_ptr, 0, point_cloud.size() - 1, point_cloud);
    Update(STATIC_ROOT_NODE);
    STATIC_ROOT_NODE->TreeSize = 0;
    Root_Node = STATIC_ROOT_NODE->left_son_ptr;
}

/**
 * @brief 执行 KD 树的最近邻搜索操作的函数
 *
 * @tparam PointType
 * @param point             搜索中心点
 * @param k_nearest         搜索最近邻点的数量
 * @param Nearest_Points    存储搜索得到的近邻点
 * @param Point_Distance    存储搜索得到的近邻点与搜索中心点的距离
 * @param max_dist          指定搜索的最远距离
 */
template <typename PointType>
void KD_TREE<PointType>::Nearest_Search(PointType point, int k_nearest, PointVector &Nearest_Points,
                                        vector<float> &Point_Distance, double max_dist)
{
    MANUAL_HEAP q(2 * k_nearest);
    q.clear();
    vector<float>().swap(Point_Distance);
    if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node)
    {
        Search(Root_Node, k_nearest, point, q, max_dist);
    }
    else
    {
        pthread_mutex_lock(&search_flag_mutex);
        while (search_mutex_counter == -1) // 重建线程正在使用搜索
        {
            pthread_mutex_unlock(&search_flag_mutex);
            usleep(1);
            pthread_mutex_lock(&search_flag_mutex);
        }
        search_mutex_counter += 1;
        pthread_mutex_unlock(&search_flag_mutex);
        Search(Root_Node, k_nearest, point, q, max_dist);
        pthread_mutex_lock(&search_flag_mutex);
        search_mutex_counter -= 1;
        pthread_mutex_unlock(&search_flag_mutex);
    }

    int k_found = min(k_nearest, int(q.size()));
    PointVector().swap(Nearest_Points);
    vector<float>().swap(Point_Distance);
    for (int i = 0; i < k_found; i++)
    {
        Nearest_Points.insert(Nearest_Points.begin(), q.top().point);
        Point_Distance.insert(Point_Distance.begin(), q.top().dist);
        q.pop();
    }
    return;
}

/**
 * @brief 在给定点的Box范围内搜索近邻点。
 *
 * @tparam PointType
 * @param Box_of_Point 给定的box范围
 * @param Storage   存储搜索结果
 */
template <typename PointType>
void KD_TREE<PointType>::Box_Search(const BoxPointType &Box_of_Point, PointVector &Storage)
{
    Storage.clear();
    Search_by_range(Root_Node, Box_of_Point, Storage);
}

/**
 * @brief 在给定点的半径范围内搜索近邻点。
 *
 * @tparam PointType
 * @param point     指定半径的圆心
 * @param radius    指定半径
 * @param Storage   存储找到的近邻点
 */
template <typename PointType>
void KD_TREE<PointType>::Radius_Search(PointType point, const float radius, PointVector &Storage)
{
    Storage.clear();
    Search_by_radius(Root_Node, point, radius, Storage);
}

/**
 * @brief 在 KD 树中添加新的数据点，如果开启了下采样，它会将新的点与下采样范围内的点进行比较，选择最近的点作为新增点。如果需要重建 KD 树，则会在操作过程中记录重建日志。
 *
 * @tparam PointType
 * @param PointToAdd    需要添加到KDTree的点向量
 * @param downsample_on 是否开启下采样，
 * @return int
 */
template <typename PointType>
int KD_TREE<PointType>::Add_Points(PointVector &PointToAdd, bool downsample_on)
{
    int NewPointSize = PointToAdd.size();
    int tree_size = size();
    BoxPointType Box_of_Point;
    PointType downsample_result, mid_point;
    bool downsample_switch = downsample_on && DOWNSAMPLE_SWITCH;
    float min_dist, tmp_dist;
    int tmp_counter = 0;
    for (int i = 0; i < PointToAdd.size(); i++)
    {
        // 开启了下采样
        if (downsample_switch)
        {
            // 计算点所在的下采样的立方体框，以确定下采样的范围。
            Box_of_Point.vertex_min[0] = floor(PointToAdd[i].x / downsample_size) * downsample_size;
            Box_of_Point.vertex_max[0] = Box_of_Point.vertex_min[0] + downsample_size;
            Box_of_Point.vertex_min[1] = floor(PointToAdd[i].y / downsample_size) * downsample_size;
            Box_of_Point.vertex_max[1] = Box_of_Point.vertex_min[1] + downsample_size;
            Box_of_Point.vertex_min[2] = floor(PointToAdd[i].z / downsample_size) * downsample_size;
            Box_of_Point.vertex_max[2] = Box_of_Point.vertex_min[2] + downsample_size;

            // 计算下采样范围的中心点
            mid_point.x = Box_of_Point.vertex_min[0] + (Box_of_Point.vertex_max[0] - Box_of_Point.vertex_min[0]) / 2.0;
            mid_point.y = Box_of_Point.vertex_min[1] + (Box_of_Point.vertex_max[1] - Box_of_Point.vertex_min[1]) / 2.0;
            mid_point.z = Box_of_Point.vertex_min[2] + (Box_of_Point.vertex_max[2] - Box_of_Point.vertex_min[2]) / 2.0;

            // 将处在下采样范围内的点添加到 Downsample_Storage 中，计算每个点到中心点的距离，并和待加入的点作比较，找到距离最近的点作为下采样的结果
            PointVector().swap(Downsample_Storage);
            Search_by_range(Root_Node, Box_of_Point, Downsample_Storage);
            min_dist = calc_dist(PointToAdd[i], mid_point);
            downsample_result = PointToAdd[i];
            for (int index = 0; index < Downsample_Storage.size(); index++)
            {
                tmp_dist = calc_dist(Downsample_Storage[index], mid_point);
                if (tmp_dist < min_dist)
                {
                    min_dist = tmp_dist;
                    downsample_result = Downsample_Storage[index];
                }
            }

            if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node)
            {
                // 如果下采样范围内有多个点，或者距离中心的最近点是待加入的点，
                if (Downsample_Storage.size() > 1 || same_point(PointToAdd[i], downsample_result))
                {
                    // 如果当前范围内已经存在点，就把这个(些)点删去
                    if (Downsample_Storage.size() > 0)
                        Delete_by_range(&Root_Node, Box_of_Point, true, true);

                    // 将待加入的点加入
                    Add_by_point(&Root_Node, downsample_result, true, Root_Node->division_axis);
                    tmp_counter++;
                }
            }
            else
            {
                // 如果正在进行第二线程的重建操作，就记录到logger中，同时以不允许重建的方式加入点
                if (Downsample_Storage.size() > 1 || same_point(PointToAdd[i], downsample_result))
                {
                    Operation_Logger_Type operation_delete, operation;
                    operation_delete.boxpoint = Box_of_Point;
                    operation_delete.op = DOWNSAMPLE_DELETE;
                    operation.point = downsample_result;
                    operation.op = ADD_POINT;
                    pthread_mutex_lock(&working_flag_mutex);
                    if (Downsample_Storage.size() > 0)
                        Delete_by_range(&Root_Node, Box_of_Point, false, true);
                    Add_by_point(&Root_Node, downsample_result, false, Root_Node->division_axis);
                    tmp_counter++;
                    if (rebuild_flag)
                    {
                        pthread_mutex_lock(&rebuild_logger_mutex_lock);
                        if (Downsample_Storage.size() > 0)
                            Rebuild_Logger.push(operation_delete);
                        Rebuild_Logger.push(operation);
                        pthread_mutex_unlock(&rebuild_logger_mutex_lock);
                    }
                    pthread_mutex_unlock(&working_flag_mutex);
                };
            }
        }
        // 没有开启下采样
        else
        {
            // 如果 Rebuild_Ptr 为空或者指向的节点不是根节点，直接插入这个点
            if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node)
            {
                Add_by_point(&Root_Node, PointToAdd[i], true, Root_Node->division_axis);
            }
            // 如果正在进行第二线程的重建操作，就将加入操作记录到logger中，并进行不允许重建的插入点操作。
            else
            {
                Operation_Logger_Type operation;
                operation.point = PointToAdd[i];
                operation.op = ADD_POINT;
                pthread_mutex_lock(&working_flag_mutex);
                Add_by_point(&Root_Node, PointToAdd[i], false, Root_Node->division_axis);
                if (rebuild_flag)
                {
                    pthread_mutex_lock(&rebuild_logger_mutex_lock);
                    Rebuild_Logger.push(operation);
                    pthread_mutex_unlock(&rebuild_logger_mutex_lock);
                }
                pthread_mutex_unlock(&working_flag_mutex);
            }
        }
    }
    return tmp_counter;
}

/**
 * @brief 将一组包围盒范围内的点（BoxPointType）重新添加到 KD 树，因此实现的事reinsertion的操作
 *
 * @tparam PointType
 * @param BoxPoints 包围盒范围内的点
 */
template <typename PointType>
void KD_TREE<PointType>::Add_Point_Boxes(vector<BoxPointType> &BoxPoints)
{
    for (int i = 0; i < BoxPoints.size(); i++)
    {
        if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node)
        {
            Add_by_range(&Root_Node, BoxPoints[i], true);
        }
        else
        {
            Operation_Logger_Type operation;
            operation.boxpoint = BoxPoints[i];
            operation.op = ADD_BOX;
            pthread_mutex_lock(&working_flag_mutex);
            Add_by_range(&Root_Node, BoxPoints[i], false);
            if (rebuild_flag)
            {
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                Rebuild_Logger.push(operation);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }
            pthread_mutex_unlock(&working_flag_mutex);
        }
    }
    return;
}

/**
 * @brief 在KDTree上删除给定点集合
 *
 * @tparam PointType
 * @param PointToDel 指定的待删除的点集合
 */
template <typename PointType>
void KD_TREE<PointType>::Delete_Points(PointVector &PointToDel)
{
    // 遍历所有点，逐个点删除
    for (int i = 0; i < PointToDel.size(); i++)
    {
        // 如果 Rebuild_Ptr 为空 或者不指向整棵树的根节点 Root_Node，直接删除要删的点
        if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node)
        {
            Delete_by_point(&Root_Node, PointToDel[i], true);
        }
        // 如果 Rebuild_Ptr 不为空且指向了 Root_Node，那么创建一个操作记录 operation，记录待删除的点以及操作类型（删除点）
        else
        {
            Operation_Logger_Type operation;
            operation.point = PointToDel[i];
            operation.op = DELETE_POINT;
            pthread_mutex_lock(&working_flag_mutex);
            Delete_by_point(&Root_Node, PointToDel[i], false);
            if (rebuild_flag)
            {
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                Rebuild_Logger.push(operation);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }
            pthread_mutex_unlock(&working_flag_mutex);
        }
    }
    return;
}

/**
 * @brief 在KDTree上删除指定包围盒范围内的点
 *
 * @tparam PointType
 * @param BoxPoints 指定要删除的box范围
 * @return int
 */
template <typename PointType>
int KD_TREE<PointType>::Delete_Point_Boxes(vector<BoxPointType> &BoxPoints)
{
    // 遍历每一个box，对每一个box单独处理
    int tmp_counter = 0;
    for (int i = 0; i < BoxPoints.size(); i++)
    {
        if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node)
        {
            tmp_counter += Delete_by_range(&Root_Node, BoxPoints[i], true, false);
        }
        else
        {
            Operation_Logger_Type operation;
            operation.boxpoint = BoxPoints[i];
            operation.op = DELETE_BOX;
            pthread_mutex_lock(&working_flag_mutex);
            tmp_counter += Delete_by_range(&Root_Node, BoxPoints[i], false, false);
            if (rebuild_flag)
            {
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                Rebuild_Logger.push(operation);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }
            pthread_mutex_unlock(&working_flag_mutex);
        }
    }
    return tmp_counter;
}

/**
 * @brief 获取需要被删除的树节点
 *
 * @tparam PointType
 * @param removed_points 返回的点向量
 */
template <typename PointType>
void KD_TREE<PointType>::acquire_removed_points(PointVector &removed_points)
{
    pthread_mutex_lock(&points_deleted_rebuild_mutex_lock);
    for (int i = 0; i < Points_deleted.size(); i++)
    {
        removed_points.push_back(Points_deleted[i]);
    }
    for (int i = 0; i < Multithread_Points_deleted.size(); i++)
    {
        removed_points.push_back(Multithread_Points_deleted[i]);
    }
    Points_deleted.clear();
    Multithread_Points_deleted.clear();
    pthread_mutex_unlock(&points_deleted_rebuild_mutex_lock);
    return;
}

/**
 * @brief 实现了KD树（K维树）数据结构中的树构建过程，通过递归沿不同轴将空间划分，构建KD树。
 *
 * @tparam PointType
 * @param root  指向当前(子)树的根节点的指针
 * @param l     表示应考虑当前(子)树中的点的范围的索引(左索引)
 * @param r     表示应考虑当前(子)树中的点的范围的索引(右索引)
 * @param Storage   包含要插入KD树的点的向量
 */
template <typename PointType>
void KD_TREE<PointType>::BuildTree(KD_TREE_NODE **root, int l, int r, PointVector &Storage)
{
    /*
    1. l > r 的出现，代表点的范围为空，因此不需要构建子树。
     */
    if (l > r)
        return;

    *root = new KD_TREE_NODE;
    InitTreeNode(*root);
    int mid = (l + r) >> 1;
    int div_axis = 0;
    int i;

    /*
    2. 根据三个维度中具有最大范围的维度来选择划分的轴。
     */
    // Find the best division Axis
    float min_value[3] = {INFINITY, INFINITY, INFINITY};
    float max_value[3] = {-INFINITY, -INFINITY, -INFINITY};
    float dim_range[3] = {0, 0, 0};
    for (i = l; i <= r; i++)
    {
        min_value[0] = min(min_value[0], Storage[i].x);
        min_value[1] = min(min_value[1], Storage[i].y);
        min_value[2] = min(min_value[2], Storage[i].z);
        max_value[0] = max(max_value[0], Storage[i].x);
        max_value[1] = max(max_value[1], Storage[i].y);
        max_value[2] = max(max_value[2], Storage[i].z);
    }
    // Select the longest dimension as division axis
    for (i = 0; i < 3; i++)
        dim_range[i] = max_value[i] - min_value[i];
    for (i = 1; i < 3; i++)
        if (dim_range[i] > dim_range[div_axis])
            div_axis = i;

    /*
    3. 递归划分
     */
    // Divide by the division axis and recursively build.
    (*root)->division_axis = div_axis;
    switch (div_axis)
    {
    case 0:
        // 将位于 Storage 容器中索引范围 [l, r] 内的元素进行部分排序，使得第 mid 个元素位于正确的位置上，
        // 其他的元素没有排序，但它的左侧元素都比它小，右侧元素都比它大
        nth_element(begin(Storage) + l, begin(Storage) + mid, begin(Storage) + r + 1, point_cmp_x);
        break;
    case 1:
        nth_element(begin(Storage) + l, begin(Storage) + mid, begin(Storage) + r + 1, point_cmp_y);
        break;
    case 2:
        nth_element(begin(Storage) + l, begin(Storage) + mid, begin(Storage) + r + 1, point_cmp_z);
        break;
    default:
        nth_element(begin(Storage) + l, begin(Storage) + mid, begin(Storage) + r + 1, point_cmp_x);
        break;
    }

    // 选择位于 mid 索引的点作为当前根的点
    (*root)->point = Storage[mid];

    // 递归构建做子树和右子树
    KD_TREE_NODE *left_son = nullptr, *right_son = nullptr;
    BuildTree(&left_son, l, mid - 1, Storage);
    BuildTree(&right_son, mid + 1, r, Storage);
    (*root)->left_son_ptr = left_son;
    (*root)->right_son_ptr = right_son;

    // 更新当前根节点信息
    Update((*root));
    return;
}

/**
 * @brief
 *
 * @tparam PointType
 * @param root  指向当前(子)树的根节点的指针
 */
template <typename PointType>
void KD_TREE<PointType>::Rebuild(KD_TREE_NODE **root)
{
    KD_TREE_NODE *father_ptr;

    // 判断当前节点的大小是否大于等于设定的 Multi_Thread_Rebuild_Point_Num，如果是，表示当前节点的大小满足多线程重建的条件。
    if ((*root)->TreeSize >= Multi_Thread_Rebuild_Point_Num)
    {
        if (!pthread_mutex_trylock(&rebuild_ptr_mutex_lock))
        {
            // 如果 Rebuild_Ptr 为空，或者当前根节点的树的节点数量大于 以Rebuild_Ptr为根节点的树的节点数量
            if (Rebuild_Ptr == nullptr || ((*root)->TreeSize > (*Rebuild_Ptr)->TreeSize))
            {
                // 将 Rebuild_Ptr 指向当前根节点
                Rebuild_Ptr = root;
            }
            pthread_mutex_unlock(&rebuild_ptr_mutex_lock);
        }
    }
    // 采用单线程重建。
    else
    {
        // 获取当前节点的父节点指针。
        father_ptr = (*root)->father_ptr;

        // 将 root 为根节点的树拉平，重建一棵KDTree，然后替换原来的以 root 为根节点的树
        int size_rec = (*root)->TreeSize;
        PCL_Storage.clear();
        flatten(*root, PCL_Storage, DELETE_POINTS_REC);
        delete_tree_nodes(root); // 删除原来的子树
        BuildTree(root, 0, PCL_Storage.size() - 1, PCL_Storage);
        if (*root != nullptr)
            (*root)->father_ptr = father_ptr;
        if (*root == Root_Node)
            STATIC_ROOT_NODE->left_son_ptr = *root;
    }
    return;
}

/**
 * @brief 在 ikd-Tree 中进行指定box范围内的节点删除
 *
 * @tparam PointType
 * @param root  指定的当前(子)树根节点
 * @param boxpoint  待删除的box范围
 * @param allow_rebuild 指定是否允许重建
 * @param is_downsample 指定是否进行下采样
 * @return int  删除的节点数量
 */
template <typename PointType>
int KD_TREE<PointType>::Delete_by_range(KD_TREE_NODE **root, BoxPointType boxpoint, bool allow_rebuild, bool is_downsample)
{
    if ((*root) == nullptr || (*root)->tree_deleted)
        return 0;
    (*root)->working_flag = true;
    Push_Down(*root);
    int tmp_counter = 0;

    /*
        进行一系列盒子范围的判断，以确定是否需要进行删除操作。
    */
    // 1. 如果盒子范围不与节点的范围相交，就返回 0 表示没有执行删除操作。
    if (boxpoint.vertex_max[0] <= (*root)->node_range_x[0] || boxpoint.vertex_min[0] > (*root)->node_range_x[1])
        return 0;
    if (boxpoint.vertex_max[1] <= (*root)->node_range_y[0] || boxpoint.vertex_min[1] > (*root)->node_range_y[1])
        return 0;
    if (boxpoint.vertex_max[2] <= (*root)->node_range_z[0] || boxpoint.vertex_min[2] > (*root)->node_range_z[1])
        return 0;

    // 2. 如果盒子范围完全包含了当前节点的范围，表示整个节点及其子树都需要被删除：
    if (boxpoint.vertex_min[0] <= (*root)->node_range_x[0] && boxpoint.vertex_max[0] > (*root)->node_range_x[1] &&
        boxpoint.vertex_min[1] <= (*root)->node_range_y[0] && boxpoint.vertex_max[1] > (*root)->node_range_y[1] &&
        boxpoint.vertex_min[2] <= (*root)->node_range_z[0] && boxpoint.vertex_max[2] > (*root)->node_range_z[1])
    {
        (*root)->tree_deleted = true;
        (*root)->point_deleted = true;
        (*root)->need_push_down_to_left = true;
        (*root)->need_push_down_to_right = true;
        tmp_counter = (*root)->TreeSize - (*root)->invalid_point_num;
        (*root)->invalid_point_num = (*root)->TreeSize;
        if (is_downsample)
        {
            (*root)->tree_downsample_deleted = true;
            (*root)->point_downsample_deleted = true;
            (*root)->down_del_num = (*root)->TreeSize;
        }
        return tmp_counter;
    }

    // 3. 盒子范围不完全包含当前节点的覆盖范围，但和当前节点的覆盖范围有交集，且包含当前节点的点，则先删除当前节点
    if (!(*root)->point_deleted &&
        boxpoint.vertex_min[0] <= (*root)->point.x && boxpoint.vertex_max[0] > (*root)->point.x &&
        boxpoint.vertex_min[1] <= (*root)->point.y && boxpoint.vertex_max[1] > (*root)->point.y &&
        boxpoint.vertex_min[2] <= (*root)->point.z && boxpoint.vertex_max[2] > (*root)->point.z)
    {
        (*root)->point_deleted = true;
        tmp_counter += 1;
        if (is_downsample)
            (*root)->point_downsample_deleted = true;
    }

    // 创建删除盒子操作的日志记录。
    Operation_Logger_Type delete_box_log;
    struct timespec Timeout;
    if (is_downsample)
        delete_box_log.op = DOWNSAMPLE_DELETE;
    else
        delete_box_log.op = DELETE_BOX;
    delete_box_log.boxpoint = boxpoint;

    // 4. 继续递归到左子树和右子树进行删除操作
    if ((Rebuild_Ptr == nullptr) || (*root)->left_son_ptr != *Rebuild_Ptr)
    {
        tmp_counter += Delete_by_range(&((*root)->left_son_ptr), boxpoint, allow_rebuild, is_downsample);
    }
    else
    {
        pthread_mutex_lock(&working_flag_mutex);
        tmp_counter += Delete_by_range(&((*root)->left_son_ptr), boxpoint, false, is_downsample);
        if (rebuild_flag)
        {
            pthread_mutex_lock(&rebuild_logger_mutex_lock);
            Rebuild_Logger.push(delete_box_log);
            pthread_mutex_unlock(&rebuild_logger_mutex_lock);
        }
        pthread_mutex_unlock(&working_flag_mutex);
    }
    if ((Rebuild_Ptr == nullptr) || (*root)->right_son_ptr != *Rebuild_Ptr)
    {
        tmp_counter += Delete_by_range(&((*root)->right_son_ptr), boxpoint, allow_rebuild, is_downsample);
    }
    else
    {
        pthread_mutex_lock(&working_flag_mutex);
        tmp_counter += Delete_by_range(&((*root)->right_son_ptr), boxpoint, false, is_downsample);
        if (rebuild_flag)
        {
            pthread_mutex_lock(&rebuild_logger_mutex_lock);
            Rebuild_Logger.push(delete_box_log);
            pthread_mutex_unlock(&rebuild_logger_mutex_lock);
        }
        pthread_mutex_unlock(&working_flag_mutex);
    }

    // 更新父节点信息
    Update(*root);

    // 判断重建相关
    if (Rebuild_Ptr != nullptr && *Rebuild_Ptr == *root && (*root)->TreeSize < Multi_Thread_Rebuild_Point_Num)
        Rebuild_Ptr = nullptr;
    bool need_rebuild = allow_rebuild & Criterion_Check((*root));
    if (need_rebuild)
        Rebuild(root);
    if ((*root) != nullptr)
        (*root)->working_flag = false;
    return tmp_counter;
}

/**
 * @brief 在KDTree上删除给定点
 *
 * @tparam PointType
 * @param root      指定(子)树根
 * @param point     给定的待删除的点
 * @param allow_rebuild 指定是否允许重建
 */
template <typename PointType>
void KD_TREE<PointType>::Delete_by_point(KD_TREE_NODE **root, PointType point, bool allow_rebuild)
{
    // 根节点为空或者整棵(子)树已经被标记为删除，直接返回
    if ((*root) == nullptr || (*root)->tree_deleted)
        return;

    (*root)->working_flag = true;
    Push_Down(*root);

    // 判断当前根节点的点是否与待删除的点匹配，如果匹配且该点没有被删除，则将该点标记为已删除，更新invalid_point_num。如果所有点都被删除，将标记该节点为已删除。
    if (same_point((*root)->point, point) && !(*root)->point_deleted)
    {
        (*root)->point_deleted = true;
        (*root)->invalid_point_num += 1;
        if ((*root)->invalid_point_num == (*root)->TreeSize)
            (*root)->tree_deleted = true;
        return;
    }

    // 如果待删除的点不是当前根节点，则继续递归到子树删除
    Operation_Logger_Type delete_log;
    struct timespec Timeout;
    delete_log.op = DELETE_POINT;
    delete_log.point = point;

    /*
    !! 新加的补救，原本有bug
     */
    if (((*root)->division_axis == 0 && point.x == (*root)->point.x) ||
        ((*root)->division_axis == 1 && point.y == (*root)->point.y) ||
        ((*root)->division_axis == 2 && point.z == (*root)->point.z))
    {
        // 如果第二线程重建的树根指针为空，或者不指向当前根节点的左孩子和右孩子
        if ((Rebuild_Ptr == nullptr) || ((*root)->left_son_ptr != *Rebuild_Ptr && (*root)->right_son_ptr != *Rebuild_Ptr))
        {
            Delete_by_point(&(*root)->left_son_ptr, point, allow_rebuild);
            Delete_by_point(&(*root)->right_son_ptr, point, allow_rebuild);
        }
        // 如果第二线程重建的树根指针不为空，并且指向当前根节点的左孩子和右孩子的其中一个
        else
        {
            // 如果第二线程重建的树根指针指向左孩子
            if ((*root)->left_son_ptr == *Rebuild_Ptr)
            {
                pthread_mutex_lock(&working_flag_mutex);
                Delete_by_point(&(*root)->left_son_ptr, point, false);
                if (rebuild_flag)
                {
                    pthread_mutex_lock(&rebuild_logger_mutex_lock);
                    Rebuild_Logger.push(delete_log);
                    pthread_mutex_unlock(&rebuild_logger_mutex_lock);
                }
                pthread_mutex_unlock(&working_flag_mutex);

                Delete_by_point(&(*root)->right_son_ptr, point, allow_rebuild);
            }
            // 如果第二线程重建的树根指针指向右孩子
            else
            {
                Delete_by_point(&(*root)->left_son_ptr, point, allow_rebuild);

                pthread_mutex_lock(&working_flag_mutex);
                Delete_by_point(&(*root)->right_son_ptr, point, false);
                if (rebuild_flag)
                {
                    pthread_mutex_lock(&rebuild_logger_mutex_lock);
                    Rebuild_Logger.push(delete_log);
                    pthread_mutex_unlock(&rebuild_logger_mutex_lock);
                }
                pthread_mutex_unlock(&working_flag_mutex);
            }
        }
    }
    else if (((*root)->division_axis == 0 && point.x < (*root)->point.x) ||
             ((*root)->division_axis == 1 && point.y < (*root)->point.y) ||
             ((*root)->division_axis == 2 && point.z < (*root)->point.z))
    {
        if ((Rebuild_Ptr == nullptr) || (*root)->left_son_ptr != *Rebuild_Ptr)
        {
            Delete_by_point(&(*root)->left_son_ptr, point, allow_rebuild);
        }
        else
        {
            pthread_mutex_lock(&working_flag_mutex);
            Delete_by_point(&(*root)->left_son_ptr, point, false);
            if (rebuild_flag)
            {
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                Rebuild_Logger.push(delete_log);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }
            pthread_mutex_unlock(&working_flag_mutex);
        }
    }
    else
    {
        if ((Rebuild_Ptr == nullptr) || (*root)->right_son_ptr != *Rebuild_Ptr)
        {
            Delete_by_point(&(*root)->right_son_ptr, point, allow_rebuild);
        }
        else
        {
            pthread_mutex_lock(&working_flag_mutex);
            Delete_by_point(&(*root)->right_son_ptr, point, false);
            if (rebuild_flag)
            {
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                Rebuild_Logger.push(delete_log);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }
            pthread_mutex_unlock(&working_flag_mutex);
        }
    }

    // 利用子节点属性更新当前根节点的属性
    Update(*root);

    // 满足条件则进行进行重建
    if (Rebuild_Ptr != nullptr && *Rebuild_Ptr == *root && (*root)->TreeSize < Multi_Thread_Rebuild_Point_Num)
        Rebuild_Ptr = nullptr;
    bool need_rebuild = allow_rebuild & Criterion_Check((*root));
    if (need_rebuild)
        Rebuild(root);

    if ((*root) != nullptr)
        (*root)->working_flag = false;
    return;
}

/**
 * @brief 将一个范围（盒子）内的点添加到 KD 树中(用于重新添加)
 *
 * @tparam PointType
 * @param root  当前(子)树根节点
 * @param boxpoint  待添加的点的范围
 * @param allow_rebuild 指定是否允许重建
 */
template <typename PointType>
void KD_TREE<PointType>::Add_by_range(KD_TREE_NODE **root, BoxPointType boxpoint, bool allow_rebuild)
{
    if ((*root) == nullptr)
        return;
    (*root)->working_flag = true;
    Push_Down(*root);

    /*
        下文中：则将当前节点的删除标志设置为下采样删除标志。
        目的是使得那些因为开启了下采样而删除的节点不会恢复成未删除状态，只恢复我们人为删除的节点。
     */
    // 检查给定的范围是否与当前节点的划分范围相交，如果不相交则返回。
    if (boxpoint.vertex_max[0] <= (*root)->node_range_x[0] || boxpoint.vertex_min[0] > (*root)->node_range_x[1])
        return;
    if (boxpoint.vertex_max[1] <= (*root)->node_range_y[0] || boxpoint.vertex_min[1] > (*root)->node_range_y[1])
        return;
    if (boxpoint.vertex_max[2] <= (*root)->node_range_z[0] || boxpoint.vertex_min[2] > (*root)->node_range_z[1])
        return;

    // 检查给定的范围是否完全包含了当前节点的划分范围，如果是，则将当前节点的删除标志设置为下采样删除标志，并将左右子节点标记为需要下推。
    if (boxpoint.vertex_min[0] <= (*root)->node_range_x[0] && boxpoint.vertex_max[0] > (*root)->node_range_x[1] &&
        boxpoint.vertex_min[1] <= (*root)->node_range_y[0] && boxpoint.vertex_max[1] > (*root)->node_range_y[1] &&
        boxpoint.vertex_min[2] <= (*root)->node_range_z[0] && boxpoint.vertex_max[2] > (*root)->node_range_z[1])
    {
        (*root)->tree_deleted = false || (*root)->tree_downsample_deleted;
        (*root)->point_deleted = false || (*root)->point_downsample_deleted;
        (*root)->need_push_down_to_left = true;
        (*root)->need_push_down_to_right = true;
        (*root)->invalid_point_num = (*root)->down_del_num;
        return;
    }

    // 当前给定范围并不完全包含整个节点的划分范围，就检查给定的范围是否包含了当前节点的点，如果是，则根据下采样删除标志设置当前节点的点删除标志。
    if (boxpoint.vertex_min[0] <= (*root)->point.x && boxpoint.vertex_max[0] > (*root)->point.x &&
        boxpoint.vertex_min[1] <= (*root)->point.y && boxpoint.vertex_max[1] > (*root)->point.y &&
        boxpoint.vertex_min[2] <= (*root)->point.z && boxpoint.vertex_max[2] > (*root)->point.z)
    {
        (*root)->point_deleted = (*root)->point_downsample_deleted;
    }

    Operation_Logger_Type add_box_log;
    struct timespec Timeout;
    add_box_log.op = ADD_BOX;
    add_box_log.boxpoint = boxpoint;

    // 递归到左子树和右子树进行re-insertion操作
    if ((Rebuild_Ptr == nullptr) || (*root)->left_son_ptr != *Rebuild_Ptr)
    {
        Add_by_range(&((*root)->left_son_ptr), boxpoint, allow_rebuild);
    }
    else
    {
        pthread_mutex_lock(&working_flag_mutex);
        Add_by_range(&((*root)->left_son_ptr), boxpoint, false);
        if (rebuild_flag)
        {
            pthread_mutex_lock(&rebuild_logger_mutex_lock);
            Rebuild_Logger.push(add_box_log);
            pthread_mutex_unlock(&rebuild_logger_mutex_lock);
        }
        pthread_mutex_unlock(&working_flag_mutex);
    }
    if ((Rebuild_Ptr == nullptr) || (*root)->right_son_ptr != *Rebuild_Ptr)
    {
        Add_by_range(&((*root)->right_son_ptr), boxpoint, allow_rebuild);
    }
    else
    {
        pthread_mutex_lock(&working_flag_mutex);
        Add_by_range(&((*root)->right_son_ptr), boxpoint, false);
        if (rebuild_flag)
        {
            pthread_mutex_lock(&rebuild_logger_mutex_lock);
            Rebuild_Logger.push(add_box_log);
            pthread_mutex_unlock(&rebuild_logger_mutex_lock);
        }
        pthread_mutex_unlock(&working_flag_mutex);
    }
    Update(*root);
    if (Rebuild_Ptr != nullptr && *Rebuild_Ptr == *root && (*root)->TreeSize < Multi_Thread_Rebuild_Point_Num)
        Rebuild_Ptr = nullptr;
    bool need_rebuild = allow_rebuild & Criterion_Check((*root));
    if (need_rebuild)
        Rebuild(root);
    if ((*root) != nullptr)
        (*root)->working_flag = false;
    return;
}

/**
 * @brief 向 KD 树中添加一个点
 *
 * @tparam PointType
 * @param root  指定的当前(子)树根节点
 * @param point 要添加的一个点
 * @param allow_rebuild 是否允许重建树
 * @param father_axis
 */
template <typename PointType>
void KD_TREE<PointType>::Add_by_point(KD_TREE_NODE **root, PointType point, bool allow_rebuild, int father_axis)
{
    // 如果当前节点为空（即到达叶子节点），则创建一个新的叶子节点，并将点添加到该节点中。根据父节点的划分轴选择下一个节点的划分轴，然后更新节点的信息并返回。
    if (*root == nullptr)
    {
        *root = new KD_TREE_NODE;
        InitTreeNode(*root);
        (*root)->point = point;
        (*root)->division_axis = (father_axis + 1) % 3;
        Update(*root);
        return;
    }

    // 如果当前节点不为空，首先设置当前节点的工作标志（working_flag）为 true，表示正在进行操作。
    (*root)->working_flag = true;

    // 创建一个操作日志对象 add_log，记录操作为添加点，并记录要添加的点。
    Operation_Logger_Type add_log;
    struct timespec Timeout;
    add_log.op = ADD_POINT;
    add_log.point = point;

    // 使用 Push_Down 函数对当前节点进行下推操作，确保当前节点及其子节点的信息是最新的。
    Push_Down(*root);

    // 根据当前节点的划分轴比较要添加的点与当前节点的划分值，决定将点添加到左子节点还是右子节点。
    if (((*root)->division_axis == 0 && point.x < (*root)->point.x) ||
        ((*root)->division_axis == 1 && point.y < (*root)->point.y) ||
        ((*root)->division_axis == 2 && point.z < (*root)->point.z))
    {
        // 如果不是正在进行第二线程重建 KD 树（即没有正在重建，或者当前节点的子节点不是正在重建的节点），则递归地在相应的子节点上进行添加操作。
        if ((Rebuild_Ptr == nullptr) || (*root)->left_son_ptr != *Rebuild_Ptr)
        {
            Add_by_point(&(*root)->left_son_ptr, point, allow_rebuild, (*root)->division_axis);
        }
        // 如果正在第二线程重建 KD 树，首先锁定 working_flag_mutex 互斥锁，然后在不允许重建的情况下添加点到相应的子节点，同时将操作日志记录到重建日志中，最后解锁互斥锁。
        else
        {
            pthread_mutex_lock(&working_flag_mutex);
            Add_by_point(&(*root)->left_son_ptr, point, false, (*root)->division_axis);
            if (rebuild_flag)
            {
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                Rebuild_Logger.push(add_log);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }
            pthread_mutex_unlock(&working_flag_mutex);
        }
    }
    else
    {
        if ((Rebuild_Ptr == nullptr) || (*root)->right_son_ptr != *Rebuild_Ptr)
        {
            Add_by_point(&(*root)->right_son_ptr, point, allow_rebuild, (*root)->division_axis);
        }
        else
        {
            pthread_mutex_lock(&working_flag_mutex);
            Add_by_point(&(*root)->right_son_ptr, point, false, (*root)->division_axis);
            if (rebuild_flag)
            {
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                Rebuild_Logger.push(add_log);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }
            pthread_mutex_unlock(&working_flag_mutex);
        }
    }

    // 更新当前节点的信息
    Update(*root);

    // 如果当前节点是正在重建的节点且节点大小小于阈值，说明不需要使用第二线程重建，将节点指针 Rebuild_Ptr 置为空
    if (Rebuild_Ptr != nullptr && *Rebuild_Ptr == *root && (*root)->TreeSize < Multi_Thread_Rebuild_Point_Num)
        Rebuild_Ptr = nullptr;

    // 根据一定的条件检查节点是否需要重建。如果需要重建，则调用 Rebuild 函数进行重建。
    bool need_rebuild = allow_rebuild & Criterion_Check((*root));
    if (need_rebuild)
        Rebuild(root);

    // 将当前节点的工作标志设为 false，表示操作完成
    if ((*root) != nullptr)
        (*root)->working_flag = false;

    return;
}

/**
 * @brief 在 k-d 树中进行最近邻搜索
 *
 * @tparam PointType
 * @param root  指定的当前(子)树根节点
 * @param k_nearest 近邻搜索的点数量
 * @param point 基准点
 * @param q 用于存储搜索结果的大根堆
 * @param max_dist 指定搜索的最远距离
 */
template <typename PointType>
void KD_TREE<PointType>::Search(KD_TREE_NODE *root, int k_nearest, PointType point, MANUAL_HEAP &q, double max_dist)
{
    if (root == nullptr || root->tree_deleted)
        return;

    // 如果 cur_dist 大于 max_dist_sqr，说明当前节点的范围盒与查询点的距离太远(整个包围盒范围内的点到基准点的距离都超过max_dist)，不需要进一步搜索，直接返回
    double cur_dist = calc_box_dist(root, point);
    double max_dist_sqr = max_dist * max_dist;
    if (cur_dist > max_dist_sqr)
        return;

    // 更新节点子节点状态信息
    int retval;
    if (root->need_push_down_to_left || root->need_push_down_to_right)
    {
        retval = pthread_mutex_trylock(&(root->push_down_mutex_lock));
        if (retval == 0)
        {
            Push_Down(root);
            pthread_mutex_unlock(&(root->push_down_mutex_lock));
        }
        else
        {
            pthread_mutex_lock(&(root->push_down_mutex_lock));
            pthread_mutex_unlock(&(root->push_down_mutex_lock));
        }
    }

    // 如果当前节点没有被删除，且查询点在当前节点范围盒内，计算点到当前节点的距离，并更新最近邻点。
    if (!root->point_deleted)
    {
        float dist = calc_dist(point, root->point);
        if (dist <= max_dist_sqr && (q.size() < k_nearest || dist < q.top().dist))
        {
            if (q.size() >= k_nearest)
                q.pop();
            PointType_CMP current_point{root->point, dist};
            q.push(current_point);
        }
    }

    int cur_search_counter;
    float dist_left_node = calc_box_dist(root->left_son_ptr, point);
    float dist_right_node = calc_box_dist(root->right_son_ptr, point);

    // 如果当前最近邻点数小于 k_nearest 或者左子节点和右子节点的距离都小于最近邻点的距离，进入递归搜索。
    if (q.size() < k_nearest || dist_left_node < q.top().dist && dist_right_node < q.top().dist)
    {
        // 根据左子节点和右子节点的距离，确定搜索顺序。如果左子节点更近，首先搜索左子节点，否则首先搜索右子节点。
        if (dist_left_node <= dist_right_node)
        {
            if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->left_son_ptr)
            {
                Search(root->left_son_ptr, k_nearest, point, q, max_dist);
            }
            // 如果 Rebuild_Ptr 指向当前子节点，首先尝试获取 search_flag_mutex
            // 如果成功，进行推入操作，然后释放锁。如果失败，表示重建线程正在进行搜索操作，等待一段时间再次尝试获取锁。
            else
            {
                pthread_mutex_lock(&search_flag_mutex);
                while (search_mutex_counter == -1) // 重建线程正在使用搜索
                {
                    pthread_mutex_unlock(&search_flag_mutex);
                    usleep(1);
                    pthread_mutex_lock(&search_flag_mutex);
                }
                search_mutex_counter += 1;
                pthread_mutex_unlock(&search_flag_mutex);
                Search(root->left_son_ptr, k_nearest, point, q, max_dist);
                pthread_mutex_lock(&search_flag_mutex);
                search_mutex_counter -= 1;
                pthread_mutex_unlock(&search_flag_mutex);
            }
            if (q.size() < k_nearest || dist_right_node < q.top().dist)
            {
                if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->right_son_ptr)
                {
                    Search(root->right_son_ptr, k_nearest, point, q, max_dist);
                }
                else
                {
                    pthread_mutex_lock(&search_flag_mutex);
                    while (search_mutex_counter == -1) // 重建线程正在使用搜索
                    {
                        pthread_mutex_unlock(&search_flag_mutex);
                        usleep(1);
                        pthread_mutex_lock(&search_flag_mutex);
                    }
                    search_mutex_counter += 1;
                    pthread_mutex_unlock(&search_flag_mutex);
                    Search(root->right_son_ptr, k_nearest, point, q, max_dist);
                    pthread_mutex_lock(&search_flag_mutex);
                    search_mutex_counter -= 1;
                    pthread_mutex_unlock(&search_flag_mutex);
                }
            }
        }
        else
        {
            if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->right_son_ptr)
            {
                Search(root->right_son_ptr, k_nearest, point, q, max_dist);
            }
            else
            {
                pthread_mutex_lock(&search_flag_mutex);
                while (search_mutex_counter == -1) // 重建线程正在使用搜索
                {
                    pthread_mutex_unlock(&search_flag_mutex);
                    usleep(1);
                    pthread_mutex_lock(&search_flag_mutex);
                }
                search_mutex_counter += 1;
                pthread_mutex_unlock(&search_flag_mutex);
                Search(root->right_son_ptr, k_nearest, point, q, max_dist);
                pthread_mutex_lock(&search_flag_mutex);
                search_mutex_counter -= 1;
                pthread_mutex_unlock(&search_flag_mutex);
            }
            if (q.size() < k_nearest || dist_left_node < q.top().dist)
            {
                if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->left_son_ptr)
                {
                    Search(root->left_son_ptr, k_nearest, point, q, max_dist);
                }
                else
                {
                    pthread_mutex_lock(&search_flag_mutex);
                    while (search_mutex_counter == -1) // 重建线程正在使用搜索
                    {
                        pthread_mutex_unlock(&search_flag_mutex);
                        usleep(1);
                        pthread_mutex_lock(&search_flag_mutex);
                    }
                    search_mutex_counter += 1;
                    pthread_mutex_unlock(&search_flag_mutex);
                    Search(root->left_son_ptr, k_nearest, point, q, max_dist);
                    pthread_mutex_lock(&search_flag_mutex);
                    search_mutex_counter -= 1;
                    pthread_mutex_unlock(&search_flag_mutex);
                }
            }
        }
    }
    // 如果当前最近邻点数 >= k_nearest 并且 左子节点和右子节点的距离不都小于最近邻点的距离
    // 那就谁小于到最近邻点的距离就到谁搜索
    else
    {
        if (dist_left_node < q.top().dist)
        {
            if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->left_son_ptr)
            {
                Search(root->left_son_ptr, k_nearest, point, q, max_dist);
            }
            else
            {
                pthread_mutex_lock(&search_flag_mutex);
                while (search_mutex_counter == -1) // 重建线程正在使用搜索
                {
                    pthread_mutex_unlock(&search_flag_mutex);
                    usleep(1);
                    pthread_mutex_lock(&search_flag_mutex);
                }
                search_mutex_counter += 1;
                pthread_mutex_unlock(&search_flag_mutex);
                Search(root->left_son_ptr, k_nearest, point, q, max_dist);
                pthread_mutex_lock(&search_flag_mutex);
                search_mutex_counter -= 1;
                pthread_mutex_unlock(&search_flag_mutex);
            }
        }
        if (dist_right_node < q.top().dist)
        {
            if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->right_son_ptr)
            {
                Search(root->right_son_ptr, k_nearest, point, q, max_dist);
            }
            else
            {
                pthread_mutex_lock(&search_flag_mutex);
                while (search_mutex_counter == -1) // 重建线程正在使用搜索
                {
                    pthread_mutex_unlock(&search_flag_mutex);
                    usleep(1);
                    pthread_mutex_lock(&search_flag_mutex);
                }
                search_mutex_counter += 1;
                pthread_mutex_unlock(&search_flag_mutex);
                Search(root->right_son_ptr, k_nearest, point, q, max_dist);
                pthread_mutex_lock(&search_flag_mutex);
                search_mutex_counter -= 1;
                pthread_mutex_unlock(&search_flag_mutex);
            }
        }
    }
    return;
}

/**
 * @brief 在给定的范围内搜索点。
 *
 * @tparam PointType
 * @param root  树的根节点
 * @param boxpoint  搜索范围
 * @param Storage   存储搜索结果
 */
template <typename PointType>
void KD_TREE<PointType>::Search_by_range(KD_TREE_NODE *root, BoxPointType boxpoint, PointVector &Storage)
{
    if (root == nullptr)
        return;
    Push_Down(root);

    // 检查给定的范围是否与当前节点的范围有重叠，如果没有重叠，则直接返回。
    if (boxpoint.vertex_max[0] <= root->node_range_x[0] || boxpoint.vertex_min[0] > root->node_range_x[1])
        return;
    if (boxpoint.vertex_max[1] <= root->node_range_y[0] || boxpoint.vertex_min[1] > root->node_range_y[1])
        return;
    if (boxpoint.vertex_max[2] <= root->node_range_z[0] || boxpoint.vertex_min[2] > root->node_range_z[1])
        return;

    // 如果给定范围完全包含了当前节点的范围，那么可以将当前节点中的点全部添加到 Storage 中
    if (boxpoint.vertex_min[0] <= root->node_range_x[0] && boxpoint.vertex_max[0] > root->node_range_x[1] &&
        boxpoint.vertex_min[1] <= root->node_range_y[0] && boxpoint.vertex_max[1] > root->node_range_y[1] &&
        boxpoint.vertex_min[2] <= root->node_range_z[0] && boxpoint.vertex_max[2] > root->node_range_z[1])
    {
        flatten(root, Storage, NOT_RECORD);
        return;
    }

    // 如果当前节点的范围与给定的范围有交集，并且当前节点中的点在给定的范围内，则将非删除状态的点添加到 Storage 中
    if (boxpoint.vertex_min[0] <= root->point.x && boxpoint.vertex_max[0] > root->point.x &&
        boxpoint.vertex_min[1] <= root->point.y && boxpoint.vertex_max[1] > root->point.y &&
        boxpoint.vertex_min[2] <= root->point.z && boxpoint.vertex_max[2] > root->point.z)
    {
        if (!root->point_deleted)
            Storage.push_back(root->point);
    }

    // 接着，检查当前节点的左子节点是否需要进行搜索，如果 Rebuild_Ptr 为空或指向的不是当前左子节点，那么直接递归调用 Search_by_range 函数。
    if ((Rebuild_Ptr == nullptr) || root->left_son_ptr != *Rebuild_Ptr)
    {
        Search_by_range(root->left_son_ptr, boxpoint, Storage);
    }
    // 如果 Rebuild_Ptr 不为空且指向的是当前左子节点，先加锁，然后递归调用 Search_by_range 函数，最后解锁。
    else
    {
        pthread_mutex_lock(&search_flag_mutex);
        Search_by_range(root->left_son_ptr, boxpoint, Storage);
        pthread_mutex_unlock(&search_flag_mutex);
    }

    // 检查当前节点的右子节点是否需要进行搜索，根据 Rebuild_Ptr 是否为空来决定直接递归调用还是先加锁再递归调用。
    if ((Rebuild_Ptr == nullptr) || root->right_son_ptr != *Rebuild_Ptr)
    {
        Search_by_range(root->right_son_ptr, boxpoint, Storage);
    }
    else
    {
        pthread_mutex_lock(&search_flag_mutex);
        Search_by_range(root->right_son_ptr, boxpoint, Storage);
        pthread_mutex_unlock(&search_flag_mutex);
    }
    return;
}

/**
 * @brief 在给定点的半径范围内搜索近邻点。
 *
 * @tparam PointType
 * @param root  指定树根
 * @param point 指定点
 * @param radius    指定半径
 * @param Storage   存储找到的近邻点
 */
template <typename PointType>
void KD_TREE<PointType>::Search_by_radius(KD_TREE_NODE *root, PointType point, float radius, PointVector &Storage)
{
    // 1. 查树节点是否为空，如果为空则直接返回。
    if (root == nullptr)
        return;

    // 2. 更新子节点的信息
    Push_Down(root);

    // 3. 计算当前节点的范围中心点
    PointType range_center;
    range_center.x = (root->node_range_x[0] + root->node_range_x[1]) * 0.5;
    range_center.y = (root->node_range_y[0] + root->node_range_y[1]) * 0.5;
    range_center.z = (root->node_range_z[0] + root->node_range_z[1]) * 0.5;

    // 4. 计算搜索点与当前节点范围中心点的欧几里得距离
    float dist = sqrt(calc_dist(range_center, point));
    // 4.1. 如果当前节点和树根的距离超出了搜索半径和树的范围半径之和，说明不会找到任何点，返回。
    if (dist > radius + sqrt(root->radius_sq))
        return;
    // 4.2 如果当前树的范围完全包含在搜索半径内，将当前树所有非标记为删除的点添加到结果中，
    if (dist <= radius - sqrt(root->radius_sq))
    {
        flatten(root, Storage, NOT_RECORD);
        return;
    }
    // 4.3 如果当前树的范围部分包含在搜索半径内，只检查当前根节点是否标记为删除，如果不是就添加到结果中
    if (!root->point_deleted && calc_dist(root->point, point) <= radius * radius)
    {
        Storage.push_back(root->point);
    }

    // 5. 对左子树进行递归搜索
    if ((Rebuild_Ptr == nullptr) || root->left_son_ptr != *Rebuild_Ptr)
    {
        Search_by_radius(root->left_son_ptr, point, radius, Storage);
    }
    else
    {
        pthread_mutex_lock(&search_flag_mutex);
        Search_by_radius(root->left_son_ptr, point, radius, Storage);
        pthread_mutex_unlock(&search_flag_mutex);
    }

    // 6. 对左子树进行递归搜索
    if ((Rebuild_Ptr == nullptr) || root->right_son_ptr != *Rebuild_Ptr)
    {
        Search_by_radius(root->right_son_ptr, point, radius, Storage);
    }
    else
    {
        pthread_mutex_lock(&search_flag_mutex);
        Search_by_radius(root->right_son_ptr, point, radius, Storage);
        pthread_mutex_unlock(&search_flag_mutex);
    }
    return;
}

/**
 * @brief 检查是否满足重建 KD 树的条件，即是否需要对 KD 树进行重建
 *
 * @tparam PointType
 * @param root  指定的(子)树节点指针
 * @return true
 * @return false
 */
template <typename PointType>
bool KD_TREE<PointType>::Criterion_Check(KD_TREE_NODE *root)
{
    // 树节点太少，无需重建
    if (root->TreeSize <= Minimal_Unbalanced_Tree_Size)
    {
        return false;
    }

    // 计算平衡度评估（balance_evaluation）和删除度评估（delete_evaluation）的值。
    // 只要有其中一个的值超过了设定阈值，就需要重建KDTree。
    float balance_evaluation = 0.0f;
    float delete_evaluation = 0.0f;
    KD_TREE_NODE *son_ptr = root->left_son_ptr;
    if (son_ptr == nullptr)
        son_ptr = root->right_son_ptr;
    delete_evaluation = float(root->invalid_point_num) / root->TreeSize;
    balance_evaluation = float(son_ptr->TreeSize) / (root->TreeSize - 1);
    if (delete_evaluation > delete_criterion_param)
    {
        return true;
    }
    if (balance_evaluation > balance_criterion_param || balance_evaluation < 1 - balance_criterion_param)
    {
        return true;
    }
    return false;
}

/**
 * @brief 将标记为下推操作的信息从父节点传递到子节点，同时更新节点的相关属性
 *
 * @tparam PointType
 * @param root  指定的树的节点
 */
template <typename PointType>
void KD_TREE<PointType>::Push_Down(KD_TREE_NODE *root)
{
    if (root == nullptr)
        return;
    Operation_Logger_Type operation;
    operation.op = PUSH_DOWN;
    operation.tree_deleted = root->tree_deleted;
    operation.tree_downsample_deleted = root->tree_downsample_deleted;
    if (root->need_push_down_to_left && root->left_son_ptr != nullptr)
    {
        // 如果当前不在重建状态或者不是在左子节点重建
        if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->left_son_ptr)
        {
            // 直接更新左子节点的删除标记和相关统计信息
            root->left_son_ptr->tree_downsample_deleted |= root->tree_downsample_deleted;
            root->left_son_ptr->point_downsample_deleted |= root->tree_downsample_deleted;
            root->left_son_ptr->tree_deleted = root->tree_deleted || root->left_son_ptr->tree_downsample_deleted;
            root->left_son_ptr->point_deleted = root->left_son_ptr->tree_deleted || root->left_son_ptr->point_downsample_deleted;
            if (root->tree_downsample_deleted)
                root->left_son_ptr->down_del_num = root->left_son_ptr->TreeSize;
            if (root->tree_deleted)
                root->left_son_ptr->invalid_point_num = root->left_son_ptr->TreeSize;
            else
                root->left_son_ptr->invalid_point_num = root->left_son_ptr->down_del_num;
            root->left_son_ptr->need_push_down_to_left = true;
            root->left_son_ptr->need_push_down_to_right = true;

            // 同时清除父节点的 need_push_down_to_left 标记
            root->need_push_down_to_left = false;
        }
        else
        {
            pthread_mutex_lock(&working_flag_mutex);
            root->left_son_ptr->tree_downsample_deleted |= root->tree_downsample_deleted;
            root->left_son_ptr->point_downsample_deleted |= root->tree_downsample_deleted;
            root->left_son_ptr->tree_deleted = root->tree_deleted || root->left_son_ptr->tree_downsample_deleted;
            root->left_son_ptr->point_deleted = root->left_son_ptr->tree_deleted || root->left_son_ptr->point_downsample_deleted;
            if (root->tree_downsample_deleted)
                root->left_son_ptr->down_del_num = root->left_son_ptr->TreeSize;
            if (root->tree_deleted)
                root->left_son_ptr->invalid_point_num = root->left_son_ptr->TreeSize;
            else
                root->left_son_ptr->invalid_point_num = root->left_son_ptr->down_del_num;
            root->left_son_ptr->need_push_down_to_left = true;
            root->left_son_ptr->need_push_down_to_right = true;

            // 如果正在进行rebuild操作，就先把操作记录到logger中
            if (rebuild_flag)
            {
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                Rebuild_Logger.push(operation);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }
            root->need_push_down_to_left = false;
            pthread_mutex_unlock(&working_flag_mutex);
        }
    }
    if (root->need_push_down_to_right && root->right_son_ptr != nullptr)
    {
        if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->right_son_ptr)
        {
            root->right_son_ptr->tree_downsample_deleted |= root->tree_downsample_deleted;
            root->right_son_ptr->point_downsample_deleted |= root->tree_downsample_deleted;
            root->right_son_ptr->tree_deleted = root->tree_deleted || root->right_son_ptr->tree_downsample_deleted;
            root->right_son_ptr->point_deleted = root->right_son_ptr->tree_deleted || root->right_son_ptr->point_downsample_deleted;
            if (root->tree_downsample_deleted)
                root->right_son_ptr->down_del_num = root->right_son_ptr->TreeSize;
            if (root->tree_deleted)
                root->right_son_ptr->invalid_point_num = root->right_son_ptr->TreeSize;
            else
                root->right_son_ptr->invalid_point_num = root->right_son_ptr->down_del_num;
            root->right_son_ptr->need_push_down_to_left = true;
            root->right_son_ptr->need_push_down_to_right = true;
            root->need_push_down_to_right = false;
        }
        else
        {
            pthread_mutex_lock(&working_flag_mutex);
            root->right_son_ptr->tree_downsample_deleted |= root->tree_downsample_deleted;
            root->right_son_ptr->point_downsample_deleted |= root->tree_downsample_deleted;
            root->right_son_ptr->tree_deleted = root->tree_deleted || root->right_son_ptr->tree_downsample_deleted;
            root->right_son_ptr->point_deleted = root->right_son_ptr->tree_deleted || root->right_son_ptr->point_downsample_deleted;
            if (root->tree_downsample_deleted)
                root->right_son_ptr->down_del_num = root->right_son_ptr->TreeSize;
            if (root->tree_deleted)
                root->right_son_ptr->invalid_point_num = root->right_son_ptr->TreeSize;
            else
                root->right_son_ptr->invalid_point_num = root->right_son_ptr->down_del_num;
            root->right_son_ptr->need_push_down_to_left = true;
            root->right_son_ptr->need_push_down_to_right = true;
            if (rebuild_flag)
            {
                pthread_mutex_lock(&rebuild_logger_mutex_lock);
                Rebuild_Logger.push(operation);
                pthread_mutex_unlock(&rebuild_logger_mutex_lock);
            }
            root->need_push_down_to_right = false;
            pthread_mutex_unlock(&working_flag_mutex);
        }
    }
    return;
}

/**
 * @brief 在构建 KD 树的过程中更新树节点的属性
 *
 * @tparam PointType
 * @param root  指定的当前(子)树根节点
 */
template <typename PointType>
void KD_TREE<PointType>::Update(KD_TREE_NODE *root)
{
    KD_TREE_NODE *left_son_ptr = root->left_son_ptr;
    KD_TREE_NODE *right_son_ptr = root->right_son_ptr;
    float tmp_range_x[2] = {INFINITY, -INFINITY};
    float tmp_range_y[2] = {INFINITY, -INFINITY};
    float tmp_range_z[2] = {INFINITY, -INFINITY};

    // Update Tree Size
    // 左右子树都不为空
    if (left_son_ptr != nullptr && right_son_ptr != nullptr)
    {
        root->TreeSize = left_son_ptr->TreeSize + right_son_ptr->TreeSize + 1;
        root->invalid_point_num = left_son_ptr->invalid_point_num + right_son_ptr->invalid_point_num + (root->point_deleted ? 1 : 0);
        root->down_del_num = left_son_ptr->down_del_num + right_son_ptr->down_del_num + (root->point_downsample_deleted ? 1 : 0);
        root->tree_downsample_deleted = left_son_ptr->tree_downsample_deleted & right_son_ptr->tree_downsample_deleted & root->point_downsample_deleted;
        root->tree_deleted = left_son_ptr->tree_deleted && right_son_ptr->tree_deleted && root->point_deleted;

        // 整棵树(包括当前根节点)都要删或者都不删
        if (root->tree_deleted || (!left_son_ptr->tree_deleted && !right_son_ptr->tree_deleted && !root->point_deleted))
        {
            tmp_range_x[0] = min(min(left_son_ptr->node_range_x[0], right_son_ptr->node_range_x[0]), root->point.x);
            tmp_range_x[1] = max(max(left_son_ptr->node_range_x[1], right_son_ptr->node_range_x[1]), root->point.x);
            tmp_range_y[0] = min(min(left_son_ptr->node_range_y[0], right_son_ptr->node_range_y[0]), root->point.y);
            tmp_range_y[1] = max(max(left_son_ptr->node_range_y[1], right_son_ptr->node_range_y[1]), root->point.y);
            tmp_range_z[0] = min(min(left_son_ptr->node_range_z[0], right_son_ptr->node_range_z[0]), root->point.z);
            tmp_range_z[1] = max(max(left_son_ptr->node_range_z[1], right_son_ptr->node_range_z[1]), root->point.z);
        }
        else
        {
            // 哪一个部分不删，那就包含哪一部分的范围
            if (!left_son_ptr->tree_deleted)
            {
                tmp_range_x[0] = min(tmp_range_x[0], left_son_ptr->node_range_x[0]);
                tmp_range_x[1] = max(tmp_range_x[1], left_son_ptr->node_range_x[1]);
                tmp_range_y[0] = min(tmp_range_y[0], left_son_ptr->node_range_y[0]);
                tmp_range_y[1] = max(tmp_range_y[1], left_son_ptr->node_range_y[1]);
                tmp_range_z[0] = min(tmp_range_z[0], left_son_ptr->node_range_z[0]);
                tmp_range_z[1] = max(tmp_range_z[1], left_son_ptr->node_range_z[1]);
            }
            if (!right_son_ptr->tree_deleted)
            {
                tmp_range_x[0] = min(tmp_range_x[0], right_son_ptr->node_range_x[0]);
                tmp_range_x[1] = max(tmp_range_x[1], right_son_ptr->node_range_x[1]);
                tmp_range_y[0] = min(tmp_range_y[0], right_son_ptr->node_range_y[0]);
                tmp_range_y[1] = max(tmp_range_y[1], right_son_ptr->node_range_y[1]);
                tmp_range_z[0] = min(tmp_range_z[0], right_son_ptr->node_range_z[0]);
                tmp_range_z[1] = max(tmp_range_z[1], right_son_ptr->node_range_z[1]);
            }
            if (!root->point_deleted)
            {
                tmp_range_x[0] = min(tmp_range_x[0], root->point.x);
                tmp_range_x[1] = max(tmp_range_x[1], root->point.x);
                tmp_range_y[0] = min(tmp_range_y[0], root->point.y);
                tmp_range_y[1] = max(tmp_range_y[1], root->point.y);
                tmp_range_z[0] = min(tmp_range_z[0], root->point.z);
                tmp_range_z[1] = max(tmp_range_z[1], root->point.z);
            }
        }
    }
    // 只有左子树不为空
    else if (left_son_ptr != nullptr)
    {
        root->TreeSize = left_son_ptr->TreeSize + 1;
        root->invalid_point_num = left_son_ptr->invalid_point_num + (root->point_deleted ? 1 : 0);
        root->down_del_num = left_son_ptr->down_del_num + (root->point_downsample_deleted ? 1 : 0);
        root->tree_downsample_deleted = left_son_ptr->tree_downsample_deleted & root->point_downsample_deleted;
        root->tree_deleted = left_son_ptr->tree_deleted && root->point_deleted;
        // 整棵树(包括当前根节点)都要删或者都不删
        if (root->tree_deleted || (!left_son_ptr->tree_deleted && !root->point_deleted))
        {
            tmp_range_x[0] = min(left_son_ptr->node_range_x[0], root->point.x);
            tmp_range_x[1] = max(left_son_ptr->node_range_x[1], root->point.x);
            tmp_range_y[0] = min(left_son_ptr->node_range_y[0], root->point.y);
            tmp_range_y[1] = max(left_son_ptr->node_range_y[1], root->point.y);
            tmp_range_z[0] = min(left_son_ptr->node_range_z[0], root->point.z);
            tmp_range_z[1] = max(left_son_ptr->node_range_z[1], root->point.z);
        }
        else
        {
            // 左子树不删，当前根节点删
            if (!left_son_ptr->tree_deleted)
            {
                tmp_range_x[0] = min(tmp_range_x[0], left_son_ptr->node_range_x[0]);
                tmp_range_x[1] = max(tmp_range_x[1], left_son_ptr->node_range_x[1]);
                tmp_range_y[0] = min(tmp_range_y[0], left_son_ptr->node_range_y[0]);
                tmp_range_y[1] = max(tmp_range_y[1], left_son_ptr->node_range_y[1]);
                tmp_range_z[0] = min(tmp_range_z[0], left_son_ptr->node_range_z[0]);
                tmp_range_z[1] = max(tmp_range_z[1], left_son_ptr->node_range_z[1]);
            }
            // 左子树删，当前根节点不删
            if (!root->point_deleted)
            {
                tmp_range_x[0] = min(tmp_range_x[0], root->point.x);
                tmp_range_x[1] = max(tmp_range_x[1], root->point.x);
                tmp_range_y[0] = min(tmp_range_y[0], root->point.y);
                tmp_range_y[1] = max(tmp_range_y[1], root->point.y);
                tmp_range_z[0] = min(tmp_range_z[0], root->point.z);
                tmp_range_z[1] = max(tmp_range_z[1], root->point.z);
            }
        }
    }
    // 只有右子树不为空
    else if (right_son_ptr != nullptr)
    {
        root->TreeSize = right_son_ptr->TreeSize + 1;
        root->invalid_point_num = right_son_ptr->invalid_point_num + (root->point_deleted ? 1 : 0);
        root->down_del_num = right_son_ptr->down_del_num + (root->point_downsample_deleted ? 1 : 0);
        root->tree_downsample_deleted = right_son_ptr->tree_downsample_deleted & root->point_downsample_deleted;
        root->tree_deleted = right_son_ptr->tree_deleted && root->point_deleted;

        // 整棵树(包括当前根节点)都要删或者都不删
        if (root->tree_deleted || (!right_son_ptr->tree_deleted && !root->point_deleted))
        {
            tmp_range_x[0] = min(right_son_ptr->node_range_x[0], root->point.x);
            tmp_range_x[1] = max(right_son_ptr->node_range_x[1], root->point.x);
            tmp_range_y[0] = min(right_son_ptr->node_range_y[0], root->point.y);
            tmp_range_y[1] = max(right_son_ptr->node_range_y[1], root->point.y);
            tmp_range_z[0] = min(right_son_ptr->node_range_z[0], root->point.z);
            tmp_range_z[1] = max(right_son_ptr->node_range_z[1], root->point.z);
        }
        else
        {
            // 右子树不删，当前根节点删
            if (!right_son_ptr->tree_deleted)
            {
                tmp_range_x[0] = min(tmp_range_x[0], right_son_ptr->node_range_x[0]);
                tmp_range_x[1] = max(tmp_range_x[1], right_son_ptr->node_range_x[1]);
                tmp_range_y[0] = min(tmp_range_y[0], right_son_ptr->node_range_y[0]);
                tmp_range_y[1] = max(tmp_range_y[1], right_son_ptr->node_range_y[1]);
                tmp_range_z[0] = min(tmp_range_z[0], right_son_ptr->node_range_z[0]);
                tmp_range_z[1] = max(tmp_range_z[1], right_son_ptr->node_range_z[1]);
            }
            // 右子树删，当前根节点不删
            if (!root->point_deleted)
            {
                tmp_range_x[0] = min(tmp_range_x[0], root->point.x);
                tmp_range_x[1] = max(tmp_range_x[1], root->point.x);
                tmp_range_y[0] = min(tmp_range_y[0], root->point.y);
                tmp_range_y[1] = max(tmp_range_y[1], root->point.y);
                tmp_range_z[0] = min(tmp_range_z[0], root->point.z);
                tmp_range_z[1] = max(tmp_range_z[1], root->point.z);
            }
        }
    }
    // 叶子结点
    else
    {
        root->TreeSize = 1;
        root->invalid_point_num = (root->point_deleted ? 1 : 0);
        root->down_del_num = (root->point_downsample_deleted ? 1 : 0);
        root->tree_downsample_deleted = root->point_downsample_deleted;
        root->tree_deleted = root->point_deleted;
        tmp_range_x[0] = root->point.x;
        tmp_range_x[1] = root->point.x;
        tmp_range_y[0] = root->point.y;
        tmp_range_y[1] = root->point.y;
        tmp_range_z[0] = root->point.z;
        tmp_range_z[1] = root->point.z;
    }

    memcpy(root->node_range_x, tmp_range_x, sizeof(tmp_range_x));
    memcpy(root->node_range_y, tmp_range_y, sizeof(tmp_range_y));
    memcpy(root->node_range_z, tmp_range_z, sizeof(tmp_range_z));

    // 计算当前节点的球形半径的平方
    float x_L = (root->node_range_x[1] - root->node_range_x[0]) * 0.5;
    float y_L = (root->node_range_y[1] - root->node_range_y[0]) * 0.5;
    float z_L = (root->node_range_z[1] - root->node_range_z[0]) * 0.5;
    root->radius_sq = x_L * x_L + y_L * y_L + z_L * z_L;

    // 更新子树的父节点
    if (left_son_ptr != nullptr)
        left_son_ptr->father_ptr = root;
    if (right_son_ptr != nullptr)
        right_son_ptr->father_ptr = root;

    // 更新平衡和删除因子
    if (root == Root_Node && root->TreeSize > 3)
    {
        KD_TREE_NODE *son_ptr = root->left_son_ptr;
        if (son_ptr == nullptr)
            son_ptr = root->right_son_ptr;
        float tmp_bal = float(son_ptr->TreeSize) / (root->TreeSize - 1);
        root->alpha_del = float(root->invalid_point_num) / root->TreeSize;
        root->alpha_bal = (tmp_bal >= 0.5 - EPSS) ? tmp_bal : 1 - tmp_bal;
    }
    return;
}

/**
 * @brief 将 KD 树节点数据(非删除状态)展平并存储到一个向量中的操作，同时还会根据不同的情况记录被删除的节点。
 *
 * @tparam PointType
 * @param root  指定的(子)树根节点
 * @param Storage 用来存储拉平树后节点的向量
 * @param storage_type  指定是否记录点，记录的点类型
 */
template <typename PointType>
void KD_TREE<PointType>::flatten(KD_TREE_NODE *root, PointVector &Storage, delete_point_storage_set storage_type)
{
    if (root == nullptr)
        return;
    Push_Down(root);

    // 收集不被标记为删除的节点
    if (!root->point_deleted)
    {
        Storage.push_back(root->point);
    }
    flatten(root->left_son_ptr, Storage, storage_type);
    flatten(root->right_son_ptr, Storage, storage_type);

    // 按需求记录要被删除的节点
    switch (storage_type)
    {
    case NOT_RECORD:
        break;
    case DELETE_POINTS_REC:
        if (root->point_deleted && !root->point_downsample_deleted)
        {
            Points_deleted.push_back(root->point);
        }
        break;
    case MULTI_THREAD_REC:
        if (root->point_deleted && !root->point_downsample_deleted)
        {
            Multithread_Points_deleted.push_back(root->point);
        }
        break;
    default:
        break;
    }
    return;
}

/**
 * @brief 递归地删除 KD 树的各个节点
 *
 * @tparam PointType
 * @param root 根指针
 */
template <typename PointType>
void KD_TREE<PointType>::delete_tree_nodes(KD_TREE_NODE **root)
{
    if (*root == nullptr)
        return;
    Push_Down(*root);
    delete_tree_nodes(&(*root)->left_son_ptr);
    delete_tree_nodes(&(*root)->right_son_ptr);

    pthread_mutex_destroy(&(*root)->push_down_mutex_lock); // 销毁当前节点的互斥锁。
    delete *root;                                          // 删除当前节点的内存，释放其占用的空间
    *root = nullptr;                                       // 将当前节点的指针设置为 nullptr，表示已删除。

    return;
}

/**
 * @brief 判断两个点是否是同一个点
 *
 * @tparam PointType
 * @param a 点1
 * @param b 点2
 * @return true
 * @return false
 */
template <typename PointType>
bool KD_TREE<PointType>::same_point(PointType a, PointType b)
{
    return (fabs(a.x - b.x) < EPSS && fabs(a.y - b.y) < EPSS && fabs(a.z - b.z) < EPSS);
}

/**
 * @brief 计算两点之间距离的平方
 *
 * @tparam PointType
 * @param a 点1
 * @param b 点2
 * @return float
 */
template <typename PointType>
float KD_TREE<PointType>::calc_dist(PointType a, PointType b)
{
    float dist = 0.0f;
    dist = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
    return dist;
}

/**
 * @brief 计算给定点到一个节点的范围盒的距离
 *
 * @tparam PointType
 * @param node  给定包围盒的中心节点
 * @param point 给定点
 * @return float 距离计算结果
 */
template <typename PointType>
float KD_TREE<PointType>::calc_box_dist(KD_TREE_NODE *node, PointType point)
{
    if (node == nullptr)
        return INFINITY;
    float min_dist = 0.0;
    if (point.x < node->node_range_x[0])
        min_dist += (point.x - node->node_range_x[0]) * (point.x - node->node_range_x[0]);
    if (point.x > node->node_range_x[1])
        min_dist += (point.x - node->node_range_x[1]) * (point.x - node->node_range_x[1]);
    if (point.y < node->node_range_y[0])
        min_dist += (point.y - node->node_range_y[0]) * (point.y - node->node_range_y[0]);
    if (point.y > node->node_range_y[1])
        min_dist += (point.y - node->node_range_y[1]) * (point.y - node->node_range_y[1]);
    if (point.z < node->node_range_z[0])
        min_dist += (point.z - node->node_range_z[0]) * (point.z - node->node_range_z[0]);
    if (point.z > node->node_range_z[1])
        min_dist += (point.z - node->node_range_z[1]) * (point.z - node->node_range_z[1]);
    return min_dist;
}

template <typename PointType>
bool KD_TREE<PointType>::point_cmp_x(PointType a, PointType b) { return a.x < b.x; }
template <typename PointType>
bool KD_TREE<PointType>::point_cmp_y(PointType a, PointType b) { return a.y < b.y; }
template <typename PointType>
bool KD_TREE<PointType>::point_cmp_z(PointType a, PointType b) { return a.z < b.z; }

// manual queue
template <typename T>
void MANUAL_Q<T>::clear()
{
    head = 0;
    tail = 0;
    counter = 0;
    is_empty = true;
    return;
}

template <typename T>
void MANUAL_Q<T>::pop()
{
    if (counter == 0)
        return;
    head++;
    head %= Q_LEN;
    counter--;
    if (counter == 0)
        is_empty = true;
    return;
}

template <typename T>
T MANUAL_Q<T>::front()
{
    return q[head];
}

template <typename T>
T MANUAL_Q<T>::back()
{
    return q[tail];
}

template <typename T>
void MANUAL_Q<T>::push(T op)
{
    q[tail] = op;
    counter++;
    if (is_empty)
        is_empty = false;
    tail++;
    tail %= Q_LEN;
}

template <typename T>
bool MANUAL_Q<T>::empty()
{
    return is_empty;
}

template <typename T>
int MANUAL_Q<T>::size()
{
    return counter;
}

template class KD_TREE<ikdTree_PointType>;
template class KD_TREE<pcl::PointXYZ>;
template class KD_TREE<pcl::PointXYZI>;
template class KD_TREE<pcl::PointXYZINormal>;