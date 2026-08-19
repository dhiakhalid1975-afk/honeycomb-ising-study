# HONEYCOMB_FINAL_HARDENING_v1_STRICT

## الغرض
حزمة post-processing داعمة نهائية للبحث:
**Critical Temperatures and Correlation-Length Exponent Identifiability in the Quenched Site-Diluted Honeycomb Ising Ferromagnet**.

هذه الحزمة **لا تشغّل Monte Carlo جديدًا، لا تغيّر v3.2.1/v3.2.1.1، ولا تستبدل القرار العلمي الأساسي**. وظيفتها عزل آخر سؤال منهجي متعلق باعتماد cubic4 على خلايا ذات ESS منخفضة خارج locked target support، ثم توليد مواد إدماج للمخطوطة من نتائج فعلية فقط.

## ما هو مقفول ولا يتغير
- البيانات الأصلية N60 وpristine read-only.
- Tc داخل كل bootstrap draw كما هو في المخرج المكتمل للحساسية غير المتناظرة.
- نفس bootstrap indices (0..1999) ونفس base seed.
- locked target support.
- symmetric feasible nu bounds.
- primary interpolation = cubic4.
- Pb residue، lattice sizes، boundary threshold 0.10، valid threshold 0.95.
- لا fallback إلى linear/PCHIP إذا فشل cubic4 بعد mask.

## الحساسية الجديدة الوحيدة
يُعاد بناء full-window corrected energy ESS من الجداول الأصلية:
`ESS_corrected = ess_energy * measurement_stride`.
الخلايا التي لا تحقق نسبة مرور 0.90 عند threshold=100 تُمنع **فقط كمصادر base interpolation**. لا تُحذف locked targets ولا يعاد بناء support بعد رؤية النتيجة.

هذه القاعدة على كامل النافذة **post-hoc sensitivity diagnostic** وليست near-Tc gate الأصلي.

## التسلسل المطلوب
1. `RUN_00_PACKAGE_TESTS.cmd`
2. `RUN_01_PRECHECK_ONLY.cmd`
3. `RUN_02_VERIFY_BASELINE_REPLAY.cmd`
4. `RUN_03_BUILD_FULL_WINDOW_ESS_AUDIT.cmd`
5. `RUN_04_RUN_ESS_DEPENDENCY_SENSITIVITY.cmd`
6. راجع `OUTPUT_FINAL_HARDENING\FINAL_ESS_DEPENDENCY_RESULT.json`
7. إذا لم يوجد HOLD: `RUN_05_BUILD_MANUSCRIPT_SUPPORT.cmd`
8. `RUN_06_PACK_OUTPUT_FOR_REVIEW.cmd`

لا تستخدم `RUN_04B_FORCE_REBUILD...` إلا إذا كنت تقصد حذف فائدة الاستئناف وإعادة الحساب النظيف؛ فهو لا يغير الأصل لكنه يعيد checkpoints الخاصة بالحساسية.

## مخرجات النشر الأساسية
- `ESS_CELL_AUDIT_REBUILT.csv`
- `LOCKED_TARGET_ESS_AUDIT.csv`
- `ESS_MASKED_DRAW_LEVEL.csv`
- `INTERPOLATION_DEPENDENCY_MAP.csv`
- `ESS_MASKED_NU_RESULTS.csv`
- `FINAL_ESS_DEPENDENCY_RESULT.json`
- `MANUSCRIPT_HARDENING_PLAN.md`
- Supplementary Table S4 (asymmetric)
- Supplementary Figure S3 (asymmetric)
- Supplementary Table S5 (ESS dependency)
- SHA-256 manifest وreview ZIP.

## قاعدة التفسير
إذا ظهرت نتيجة قد تغيّر القرار (القناتان تمران البوابات وتستبعدان nu=1 على الجانب نفسه) فالحزمة لا "تصلح" النتيجة؛ بل تضع `POTENTIAL_DECISION_IMPACT_REQUIRES_FULL_REVIEW` ويجب إيقاف إدماج المخطوطة حتى المراجعة.

## ملاحظة عن asymmetric audit
المخرج المكتمل مضمّن كمرجع immutable. العدد الصحيح للمخفف هو **5/6** upper endpoints عند 1.45، والاستثناء p=0.90 Binder (upper=1.273515457367...). هذه الحساسية لا تستبدل symmetric primary analysis.

## حماية سلامة الحزمة نفسها
قبل قراءة بيانات المشروع، يفحص `RUN_01` بصمات SHA-256 للحمولة البرمجية والمرجعية المقفلة عبر `IMMUTABLE_PACKAGE_SHA256_MANIFEST.csv`. الملفات المحلية القابلة للتعديل (`USER_CONFIG.json` ومسار Python/المشروع) مستثناة عمدًا من هذا القفل. أي تغيير غير مصرح في الكود أو المرجع يؤدي إلى FAIL-CLOSED.

## تقرير المراجعة
بعد `RUN_04` يُنشأ `FINAL_ESS_DEPENDENCY_REPORT.md` إلى جانب JSON والجداول التفصيلية. لا تعتمد أي صياغة للمخطوطة قبل قراءة `FINAL_ESS_DEPENDENCY_RESULT.json` و`MANUSCRIPT_INTEGRATION_GUARD.json`.

## نطاق إغلاق ESS
قرار `strong_ess_dependency_closure_supported` يخص الحالات المخففة الثلاث فقط. يتم تشغيل pristine أيضًا كـpost-hoc stress diagnostic للشفافية، لكنه لا يعيد تعريف معايرة p=1 الأصلية لأن full-window 100/90% rule ليست بوابة الدراسة الأصلية. يظهر هذا بوضوح في التقرير الكامل، بينما جدول المخطوطة S5 يقتصر على الحالات المخففة.
