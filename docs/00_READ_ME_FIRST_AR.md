# FGT_PUBLICATION_FREEZE_v3_2_1_1

هذه لقطة التجميد العلمي المنظمة للعمل النهائي.

## قاعدة الاعتماد
- الاعتماد العلمي على الملفات الفردية المقفلة + SHA-256 + provenance.
- الحزمة الأصلية v3.2.1 وSPEC_LOCK محفوظان دون تعديل في `01_ORIGINAL_IMMUTABLE`.
- تصحيح طبقة القرار v3.2.1.1 منفصل في `02_FINAL_DECISION_CORRECTION` ولا يعيد حساب Monte Carlo أو bootstrap أو Tc أو nu.
- أدلة ESS/time-series/equilibration منفصلة في `03_FINAL_AUDIT_EVIDENCE`.
- كود الأشكال ومجموعة الأشكال المضمنة فعلياً في المخطوطة النهائية موجودة في `04_FINAL_FIGURE_CODE_AND_FIGURES`.
- المخطوطة النهائية وهيكلها ومراجعها في `05_PAPER_STRUCTURE_AND_REFERENCES`.
- الحزم الست التاريخية بقيت غير معدلة في مكتبة المشروع، وجرى تثبيت أسماءها ومساراتها ومعرفاتها وأحجامها في `90_HISTORICAL_ARCHIVE/HISTORICAL_ARCHIVE_POINTERS.csv`. لم يُنشأ أي ZIP بديل مزيف عندما تعذر تصدير raw bytes.

## منع الاستبدال
- لا يستبدل `TABLE_FINAL_DECISIONS_v321_1.csv` بجدول v321 الأقدم.
- لا يستخدم legacy fixed-nu susceptibility-shift Tc كـTc inferential نهائي أو لاختبار nu.
- لا تتحول robustness envelope إلى confidence interval.
- لا تصدر دعوى universality جديدة من الحالات المخففة لأن nu غير قابل للتحديد ضمن L=40–120.

راجع `FREEZE_SHA256_MANIFEST.csv` و`FREEZE_INTEGRITY_REPORT.json` للتحقق الآلي.
