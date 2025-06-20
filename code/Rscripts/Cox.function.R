###################################################################################
#函数名称：Cox.function
#' @description 对临床变量进行单多因素COX分析，以表格形式输出结果
#' @param time  数值型向量，代表样本生存事件对应的时间
#' @param event 数值型向量，代表样本生存事件对应的结局，0 = alive，1 = dead
#' @param clinical.data  样本临床数据信息，至少包含两个内容：（1）样本ID；（2）样本的临床变量信息（如：年龄，tumor stage）
#' @param clinical.variate 数值型向量，用来指明clinical.data中需要分析的临床变量所在列。默认将从第四列开始的变量作为临床变量
#' @return 返回一个数据框，包含单多因素COX分析的结果信息,信息包括单多因素的HR（95% CI for HR），p value
Cox.function <- function(time, event, clinical.data, clinical.variate = NULL)
{
  ###设置工作环境
  options(stringsAsFactors = FALSE, warn = -1);
  suppressPackageStartupMessages(require(survival));

  ###判断协变量类型：数值型（num.covariate），非数值型(chara.covariate)。便于后续输出
  if(is.null(clinical.variate))
  {
      covariates  <- colnames(clinical.data)[-c(1:3)];
  }
  if(is.numeric(clinical.variate))
  {
      covariates  <- colnames(clinical.data)[clinical.variate];
  }
  num.variate <- NULL
  for(i in covariates)
  {
      if(is.numeric(clinical.data[, i]))
      {
          num.variate <- append(num.variate, i)
      }
  }
  chara.variate <- setdiff(covariates, num.variate)

  ####单因素cox分析函数：univariate.cox
  univariate.cox <- function(data, num, chara)
  {
    #STEP1:构建单因素分析的对象
    univ_formulas <- sapply(covariates,
                            function(x) as.formula(paste('Surv(time, event)~', x)));
    
    #STEP2:单因素Cox分析
    univ_models <- lapply(univ_formulas, function(x){coxph(x, data = data)});

    #STEP3:提取有用信息
    univ_results <- lapply(univ_models, function(x)
                            {                             
                             tmp <-summary(x);
                             
                             #提取p值，保留两位有效数字
                             p.value <- round(tmp$coefficients[ ,5], digits = 4);
                             p.value[which(p.value < 0.0001)] <- "<0.0001";
                             
                             #提取beta值，这里的coefficients为矩阵，但是只有一行，所以可以这样取值
                             #beta <- round(tmp$coefficients[ ,1], digits = 4);
                             
                             #提取风险比
                             HR <- round(tmp$coefficients[ ,2], digits = 4);
                             
                             #提取95%置信区间上下界
                             HR.confint.lower <- round(tmp$conf.int[,"lower .95"], 4);
                             HR.confint.upper <- round(tmp$conf.int[,"upper .95"], 4);    

                             #合并风险比HR和置信区间为一个内容
                             HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")");
                             
                             variate <- rownames(tmp$coefficients);
                             
                             #将所有值合并在一个矩阵中
                             all.data <- as.data.frame(cbind(variate, HR, p.value));
                           }
                          );
    
    #STEP4:标准化输出格式
    for(i in num)
    {
        tmp <- univ_results[[i]];
        tmp$type <- " ";
        univ_results[[i]] <- tmp; 
    }

    for(i in chara)
    {
        tmp <- univ_results[[i]]
        tmp.variate <- substr(tmp$variate, start = 0, stop = nchar(i));
        tmp.type <- sub(pattern = i, replacement = "", tmp$variate);
        tmp$variate <- tmp.variate;
        tmp$type <- paste(tmp.type, "VS", setdiff(unique(data[, i]), tmp.type));
        univ_results[[i]] <- tmp; 
    }

    #STEP5: list转化为数据框输出 	
    univ.result <- do.call(rbind.data.frame, univ_results);
    univ.result <- univ.result[,c(1,4,2,3)];
    colnames(univ.result) <- c('variate', 'type', 'univ HR (95% CI for HR)', 'univ p value')
    rownames(univ.result) <- NULL;
    return(univ.result)  
  };

  ####多因素cox分析函数：multivariate.cox
  multivariate.cox <- function(data, num, chara)
  {
    options(stringsAsFactors = FALSE);

    #STEP1:直接对所有选中的协变量进行多因素分析
    multiv_formula <- as.formula(paste("Surv(time, event)~", paste(covariates, collapse="+"), sep=""));
    multiv_model    <- coxph(multiv_formula, data = data);
    
    #STEP2:提取有用信息
    tmp <- summary(multiv_model);
    
    #提取p值，保留两位有效数字
    p.value <- round(tmp$coefficients[ ,5], digits = 4);
    p.value[which(p.value < 0.0001)] <- "<0.0001"; 
    
    #提取beta值，这里得到的coefficients为矩阵，且有多行（每行对应一个协变量）
    #beta <- round(tmp$coefficients[ ,1], digits = 4);
    
    #提取风险比
    HR <- round(tmp$coefficients[ ,2], digits = 4);
    
    #提取95%置信区间上下界
    HR.confint.lower <- round(tmp$conf.int[ ,"lower .95"], digits = 4);
    HR.confint.upper <- round(tmp$conf.int[ ,"upper .95"],digits = 4);
    
    #合并风险比HR和置信区间为一个内容
    HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")");
    
    variate <- rownames(tmp$coefficients);
    
    #整合输出内容为data.frame
    multiv_result <- as.data.frame(cbind(variate, HR, p.value));

    #STEP3:新建数据框储存多因素结果
    multiv.result <- NULL;

    for(i in num)
    {
        n.row <- grep(pattern = i, multiv_result$variate);
        tmp <- multiv_result[n.row, ];
        tmp$type <- " ";
        multiv.result <- rbind(multiv.result,tmp);
    }

    for(i in chara)
    {
        n.row <- grep(pattern = i, multiv_result$variate);
        tmp <- multiv_result[n.row, ];
        tmp.variate <- substr(tmp$variate, start = 0, stop = nchar(i));
        tmp.type <- sub(pattern = i, replacement = "", tmp$variate);
        tmp$variate <- tmp.variate;
        tmp$type <- paste(tmp.type, "VS", setdiff(unique(data[, i]), tmp.type));
        multiv.result <- rbind(multiv.result,tmp);         
    }
    multiv.result <- multiv.result[,c(1,4,2,3)]
    colnames(multiv.result) <- c("variate", "type","multiv HR (95% CI for HR)", "multiv p value");
    rownames(multiv.result) <- NULL;

    return(multiv.result);
  };

  ###运行上述函数
  UniCoxPH <- univariate.cox(data = clinical.data, num = num.variate, chara = chara.variate);
  MultiCoxPH <- multivariate.cox(data = clinical.data, num = num.variate, chara = chara.variate);

  ###合并两个数据框
  cox.result <- merge(UniCoxPH, MultiCoxPH, by = c("variate", "type"), all = T);
  colnames(cox.result) <- c("variate", " ", "univ HR (95% CI for HR)", "univ p value","multiv HR (95% CI for HR)", "multiv p value");

  ###更改表格格式
  for(i in chara.variate)
  {
  tmp.row <- which(cox.result[,1] == i);
  tmp.vec <- c(i, rep(" ", times = 5));
  cox.result <- rbind(cox.result[1 : (tmp.row-1),], tmp.vec, cox.result[tmp.row : nrow(cox.result),]); 
  };
  cox.result[duplicated(cox.result[,1]),1] <- " ";
  rownames(cox.result) <- 1:nrow(cox.result);

  return(cox.result)
}


