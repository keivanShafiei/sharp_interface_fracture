#ifndef FIELD_TRANSFER_HPP
#define FIELD_TRANSFER_HPP

#include "mesh.hpp"
#include <vector>

/**
 * @brief کلاس انتقال میدان برای نگاشت نتایج بین مش‌های مختلف
 * * این کلاس وظیفه دارد مقادیر جابجایی (Displacement) را از یک مش قدیمی (قبل از بازسازی)
 * به مش جدید (بعد از بازسازی در اطراف ترک) منتقل کند تا پایداری حلگر حفظ شود.
 */
class FieldTransfer {
public:
    FieldTransfer() = default;

    /**
     * @brief انتقال بردار جابجایی از مش قدیمی به جدید
     * * @param old_mesh مش قبلی که حل روی آن انجام شده
     * @param new_mesh مش جدید ایجاد شده توسط Remesher
     * @param u_old بردار جابجایی مش قبلی
     * @param u_new بردار جابجایی مش جدید (خروجی)
     * @return true در صورت موفقیت‌آمیز بودن عملیات
     */
    bool transfer_displacement(const Mesh& old_mesh,
                               const Mesh& new_mesh,
                               const std::vector<double>& u_old,
                               std::vector<double>& u_new) const;

private:
    /**
     * @brief تابع کمکی برای محاسبه مختصات باریسنتریک
     * برای پیدا کردن موقعیت دقیق گره جدید در المان قدیمی
     */
    bool compute_barycentric(double px, double py,
                             double x1, double y1, 
                             double x2, double y2, 
                             double x3, double y3,
                             double& w1, double& w2, double& w3) const;
};

#endif // FIELD_TRANSFER_HPP
