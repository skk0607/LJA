//
// Created by anton on 19.12.2019.
//

#include "contigs.hpp"

bool StringContig::homopolymer_compressing = false;
size_t StringContig::min_dimer_to_compress = 1000000000;
size_t StringContig::max_dimer_size = 1000000000;
size_t StringContig::dimer_step = 1;

/**
 * @brief 压缩字符串连续序列
 *
 * 如果未启用同聚物压缩，则不进行任何操作。
 * 否则，删除序列中的重复元素，并进行同聚物压缩。
 *
 * @param void
 *
 * @return void
 */
void StringContig::compress() {
    if(!homopolymer_compressing)
        return;

    // 如果homopolymer_compressing为false，则直接返回，不进行压缩操作
    //主要目的是从 std::vector 中删除所有重复的连续元素。
    seq.erase(std::unique(seq.begin(), seq.end()), seq.end());

    // 使用std::unique函数去除seq中连续的重复元素
    VERIFY(min_dimer_to_compress <= max_dimer_size);
    VERIFY(min_dimer_to_compress >= 4);
    VERIFY(dimer_step == 1);

    // 验证min_dimer_to_compress是否在max_dimer_size的范围内
    // 验证min_dimer_to_compress是否大于等于4
    // 验证dimer_step是否等于1
    if(min_dimer_to_compress >= seq.size())
        return;

    // 如果min_dimer_to_compress大于等于seq的大小，则直接返回，不进行压缩操作
    size_t cur = 2;
    size_t at_len = 2;

    // 定义当前位置和连续相同元素的长度
    for(size_t i = 2; i <= seq.size(); i++) {
        if (i < seq.size() && seq[i] == seq[cur - 2]) {
            seq[cur] = seq[i];
            cur++;
            at_len += 1;
        } else {
            // 如果当前元素与前一个元素相同，则更新当前位置和连续相同元素的长度
            if(at_len > min_dimer_to_compress) {
                // 如果连续相同元素的长度大于min_dimer_to_compress
                if(at_len > max_dimer_size) {
                    // 如果连续相同元素的长度大于max_dimer_size
                    cur -= (at_len - max_dimer_size) / 2 * 2;
                    at_len -= (at_len - max_dimer_size) / 2 * 2;
                }
                size_t corr_at_len = at_len / dimer_step * dimer_step;
                VERIFY(corr_at_len == at_len);
                cur -= (at_len - corr_at_len) / 2 * 2;
                at_len -= (at_len - corr_at_len) / 2 * 2;
            }
            at_len = 2;

            // 调整当前位置和连续相同元素的长度
            if(i < seq.size()) {
                seq[cur] = seq[i];
                cur++;
            }
        }
    }
    seq.erase(seq.begin() + cur, seq.end());

    // 移除seq中从当前位置到末尾的元素
}
