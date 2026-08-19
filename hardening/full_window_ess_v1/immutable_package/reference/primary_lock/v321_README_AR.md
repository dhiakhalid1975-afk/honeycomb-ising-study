# FGT Correction-Aware Critical Audit v3.2.1

## الغرض
هذه حزمة **Post-processing = معالجة لاحقة** لبيانات N60 المقبولة في مشروع:
`FGT_Dilution_Study_Code_v2.4.0_PUBLICATION_STRICT`.

لا تعيد Monte Carlo، ولا تغيّر ملفات المشروع الأصلية، ولا تحقن قيم DFT. هدفها معالجة المشكلة المنهجية التي ظهرت في التدقيق السابق: فصل تقدير `Tc` عن فرض `nu=1`، ثم اختبار التحجيم الحرج مع تشخيص صريح لتصحيحات الحجم المحدود.

## ماذا تفعل علمياً؟
- تدقق قلب Monte Carlo الأصلي ديناميكياً وتثبت SHA-256 لملفاته الحرجة.
- تستخرج `Tc` من تقاطعات Binder و `xi/L` بدل تقييده بفاصل Tc الناتج من fixed-nu shift.
- تلائم `nu` ببحث أحادي البعد بعد تثبيت Tc المستقل.
- تعيد تقدير Tc ثم nu داخل كل سحبة quenched bootstrap.
- تحسب paired `L_min` drift بين `L=40..120` و `L=60..120`.
- تميز بين bootstrap CI وبين robustness envelope المنهجي.
- تستخدم p=1 كمعايرة خارجية ضد Tc التحليلي ومرجعيات Ising للأسس.
- تمنع تلقائياً الادعاء الزائد عند ظهور drift أو boundary/validity problems أو فشل المعايرة.

## الترتيب الصحيح للتشغيل على Windows
ضع الحزمة في مجلد مستقل ولا تخلطها مع v3.1.2. شغّل من داخل هذا المجلد:

```text
00_INSTALL_AND_VALIDATE_ENV.cmd
01_SYNTHETIC_CHALLENGE_4_WORKERS.cmd
02_BOOTSTRAP_CONVERGENCE_4_WORKERS.cmd
03_RUN_OR_RESUME_REAL_AUDIT_4_WORKERS.cmd
06_VALIDATE_RELEASE.cmd
```

`04_STATUS.cmd` يعرض الحالة فقط، و`05_REBUILD_FIGURES_ONLY.cmd` يعيد الرسوم دون إعادة bootstrap.

عند أول تشغيل لـ`00` ألصق المسار الكامل لجذر المشروع الأصلي. بعد ذلك يُحفظ في `PROJECT_ROOT.local.txt` داخل نسخة الحزمة نفسها.

## أربعة عمال والاستئناف
الافتراضي أربعة workers. على Windows يستخدم backend خيطي لتجنب مشاكل console handles. بذور bootstrap تعتمد على case/bootstrap index ولا تعتمد على ترتيب العامل.

إذا انطفأ الجهاز أثناء `02` أو `03`، شغّل نفس ملف CMD مرة أخرى. كل chunk مكتمل لا يُقبل إلا بعد تطابق run signature وSHA-256.

## مساحة القرص
على Windows تحفظ checkpoints في مسار قصير تحت `%USERPROFILE%\FGTCSA` افتراضياً. قبل bootstrap يوجد I/O probe وفحص مساحة حرة. عند نقص المساحة يمكن تحديد قرص آخر قبل التشغيل:

```cmd
set FGT_CSAUDIT_WORK_ROOT=D:\FGTCSA_WORK
```

ثم أعد نفس الخطوة. لا تحذف checkpoints الصحيحة إذا كنت تريد الاستئناف.

## المخرجات النهائية
تكتب فقط إلى:
`results\publication_strict_phase3\final_csaudit_v3_2_1_correction_aware\`

أهم الجداول:
- `TABLE_CENTRAL_RG_TC_AND_NU.csv`
- `TABLE_SUPPORT_BALANCE_AUDIT.csv`
- `TABLE_TC_RG_CROSSING_BOOTSTRAP.csv`
- `TABLE_NU_CORRECTION_AWARE.csv`
- `TABLE_PAIRED_LMIN_DRIFT.csv`
- `TABLE_EXPONENT_RATIO_CONSISTENCY.csv`
- `TABLE_CHANNEL_DELTA_NU.csv`
- `TABLE_FINAL_DECISIONS_v321.csv`
- `CLAIM_SCOPE_TABLE_v321.csv`

الرسوم تحفظ PDF/SVG وPNG بدقة 600 dpi.

## لغة النتائج المسموحة
القرارات الممكنة:
- `NU1_COMPATIBLE`
- `INCONCLUSIVE_CORRECTION_DOMINATED`
- `INCONCLUSIVE_LIMITED_RANGE`
- `EVIDENCE_AGAINST_NU1`

حتى القرار الأخير لا يساوي إثبات universality class جديدة.

## مهم
نجاح اختبارات البرمجة لا يعني أن نتيجة nu ستكون حاسمة. مع L=40..120 قد تكون النتيجة العلمية الصحيحة هي `INCONCLUSIVE` بسبب تصحيحات الحجم المحدود/اللوغاريتمية. هذا سلوك مقصود fail-closed وليس فشلاً برمجياً.
