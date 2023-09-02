# This part of the code plot both return series and volatility together

# US Spillover to other countries first set
US_to_Other = theta1mat[,1]
state = data.frame(US_to_Other)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover from US to Other countries") +ylim(0,1.0)

# US Spillover from other countries second set

US_from_Other = theta1mat[,2]
state = data.frame(US_from_Other)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover to US from Other countries") +ylim(0,1.75)

# US Spillover to other countries: good volatility # third

US_to_Other_gd = theta1mat[,3]
US_to_Other_gd0 = theta1mat[,5]
state = data.frame(US_to_Other_gd,US_to_Other_gd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover from US to Other countries")+ylim(0,1.0)

# US Spillover from other countries: good volatility fourth

US_from_Other_gd = theta1mat[,4]
US_from_Other_gd0 = theta1mat[,6]
state = data.frame(US_from_Other_gd,US_from_Other_gd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover to US from Other countries")+ylim(0,1.75)

# US Spillover to other countries: bad volatility fifth

US_to_Other_bd = theta1mat[,7]
US_to_Other_bd0 = theta1mat[,9]
state = data.frame(US_to_Other_bd,US_to_Other_bd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover from US to Other countries")+ylim(0,1.0)

# US Spillover from other countries: bad volatility sixth

US_from_Other_bd = theta1mat[,8]
US_from_Other_bd0 = theta1mat[,10]
state = data.frame(US_from_Other_bd,US_from_Other_bd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover to US from Other countries")+ylim(0,1.75)

# Canada

# This part of the code plot both return series and volatility together

# Canada Spillover to other countries first set
Canada_to_Other = theta2mat[,1]
state = data.frame(Canada_to_Other)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover from Canada to Other countries")+ylim(0,1.0)

# Canada Spillover from other countries second set

Canada_from_Other = theta2mat[,2]
state = data.frame(Canada_from_Other)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover to Canada from Other countries")+ylim(0,1.75)

# Canada Spillover to other countries: good volatility # third

Canada_to_Other_gd = theta2mat[,3]
Canada_to_Other_gd0 = theta2mat[,5]
state = data.frame(Canada_to_Other_gd,Canada_to_Other_gd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Good Volatility Spillover from Canada to Other countries")+ylim(0,1.0)

# Canada Spillover from other countries: good volatility fourth

Canada_from_Other_gd = theta2mat[,4]
Canada_from_Other_gd0 = theta2mat[,6]
state = data.frame(Canada_from_Other_gd,Canada_from_Other_gd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Good Volatility Spillover to Canada from Other countries")+ylim(0,1.75)

# Canada Spillover to other countries: bad volatility fifth

Canada_to_Other_bd = theta2mat[,7]
Canada_to_Other_bd0 = theta2mat[,9]
state = data.frame(Canada_to_Other_bd,Canada_to_Other_bd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Bad Volatility Spillover from Canada to Other countries")+ylim(0,1.0)

# Canada Spillover from other countries: bad volatility sixth

Canada_from_Other_bd = theta2mat[,8]
Canada_from_Other_bd0 = theta2mat[,10]
state = data.frame(Canada_from_Other_bd,Canada_from_Other_bd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Bad Volatility Spillover to Canada from Other countries")+ylim(0,1.75)

# Japan

# This part of the code plot both return series and volatility together

# Japan Spillover to other countries first set
Japan_to_Other = theta3mat[,1]
state = data.frame(Japan_to_Other)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover from Japan to Other countries")+ylim(0,1.0)

# Japan Spillover from other countries second set

Japan_from_Other = theta3mat[,2]
state = data.frame(Japan_from_Other)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover to Japan from Other countries")+ylim(0,1.75)

# Japan Spillover to other countries: good volatility # third

Japan_to_Other_gd = theta3mat[,3]
Japan_to_Other_gd0 = theta3mat[,5]
state = data.frame(Japan_to_Other_gd,Japan_to_Other_gd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Good Volatility Spillover from Japan to Other countries")+ylim(0,1.0)

# Japan Spillover from other countries: good volatility fourth

Japan_from_Other_gd = theta3mat[,4]
Japan_from_Other_gd0 = theta3mat[,6]
state = data.frame(Japan_from_Other_gd,Japan_from_Other_gd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Good Volatility Spillover to Japan from Other countries")+ylim(0,1.75)

# Japan Spillover to other countries: bad volatility fifth

Japan_to_Other_bd = theta3mat[,7]
Japan_to_Other_bd0 = theta3mat[,9]
state = data.frame(Japan_to_Other_bd,Japan_to_Other_bd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Bad Volatility Spillover from Japan to Other countries")+ylim(0,1.0)

# Japan Spillover from other countries: bad volatility sixth

Japan_from_Other_bd = theta3mat[,8]
Japan_from_Other_bd0 = theta3mat[,10]
state = data.frame(Japan_from_Other_bd,Japan_from_Other_bd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Bad Volatility Spillover to Japan from Other countries")+ylim(0,1.75)

# Italy

# This part of the code plot both return series and volatility together

# Italy Spillover to other countries first set
Italy_to_Other = theta4mat[,1]
state = data.frame(Italy_to_Other)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover from Italy to Other countries")+ylim(0,1.0)

# Italy Spillover from other countries second set

Italy_from_Other = theta4mat[,2]
state = data.frame(Italy_from_Other)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover to Italy from Other countries")+ylim(0,1.75)

# Italy Spillover to other countries: good volatility # third

Italy_to_Other_gd = theta4mat[,3]
Italy_to_Other_gd0 = theta4mat[,5]
state = data.frame(Italy_to_Other_gd,Italy_to_Other_gd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Good Volatility Spillover from Italy to Other countries")+ylim(0,1.0)

# Italy Spillover from other countries: good volatility fourth

Italy_from_Other_gd = theta4mat[,4]
Italy_from_Other_gd0 = theta4mat[,6]
state = data.frame(Italy_from_Other_gd,Italy_from_Other_gd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Good Volatility Spillover to Italy from Other countries")+ylim(0,1.75)

# Italy Spillover to other countries: bad volatility fifth

Italy_to_Other_bd = theta4mat[,7]
Italy_to_Other_bd0 = theta4mat[,9]
state = data.frame(Italy_to_Other_bd,Italy_to_Other_bd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Bad Volatility Spillover from Italy to Other countries")+ylim(0,1.0)

# Italy Spillover from other countries: bad volatility sixth

Italy_from_Other_bd = theta4mat[,8]
Italy_from_Other_bd0 = theta4mat[,10]
state = data.frame(Italy_from_Other_bd,Italy_from_Other_bd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Bad Volatility Spillover to Italy from Other countries")+ylim(0,1.75)


#Germany
# This part of the code plot both return series and volatility together

# Germany Spillover to other countries first set
Germany_to_Other = theta5mat[,1]
state = data.frame(Germany_to_Other)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover from Germany to Other countries")+ylim(0,1.0)

# Germany Spillover from other countries second set

Germany_from_Other = theta5mat[,2]
state = data.frame(Germany_from_Other)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover to Germany from Other countries")+ylim(0,1.75)

# Germany Spillover to other countries: good volatility # third

Germany_to_Other_gd = theta5mat[,3]
Germany_to_Other_gd0 = theta5mat[,5]
state = data.frame(Germany_to_Other_gd,Germany_to_Other_gd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Good Volatility Spillover from Germany to Other countries")+ylim(0,1.0)

# Germany Spillover from other countries: good volatility fourth

Germany_from_Other_gd = theta5mat[,4]
Germany_from_Other_gd0 = theta5mat[,6]
state = data.frame(Germany_from_Other_gd,Germany_from_Other_gd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Good Volatility Spillover to Germany from Other countries")+ylim(0,1.75)

# Germany Spillover to other countries: bad volatility fifth

Germany_to_Other_bd = theta5mat[,7]
Germany_to_Other_bd0 = theta5mat[,9]
state = data.frame(Germany_to_Other_bd,Germany_to_Other_bd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Bad Volatility Spillover from Germany to Other countries")+ylim(0,1.0)

# Germany Spillover from other countries: bad volatility sixth

Germany_from_Other_bd = theta5mat[,8]
Germany_from_Other_bd0 = theta5mat[,10]
state = data.frame(Germany_from_Other_bd,Germany_from_Other_bd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Bad Volatility Spillover to Germany from Other countries")+ylim(0,1.75)


#France
# This part of the code plot both return series and volatility together

# France Spillover to other countries first set
France_to_Other = theta6mat[,1]
state = data.frame(France_to_Other)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover from France to Other countries")+ylim(0,1.0)

# France Spillover from other countries second set

France_from_Other = theta6mat[,2]
state = data.frame(France_from_Other)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Volatility Spillover to France from Other countries")+ylim(0,1.75)

# France Spillover to other countries: good volatility # third

France_to_Other_gd = theta6mat[,3]
France_to_Other_gd0 = theta6mat[,5]
state = data.frame(France_to_Other_gd,France_to_Other_gd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Good Volatility Spillover from France to Other countries")+ylim(0,1.0)

# France Spillover from other countries: good volatility fourth

France_from_Other_gd = theta6mat[,4]
France_from_Other_gd0 = theta6mat[,6]
state = data.frame(France_from_Other_gd,France_from_Other_gd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Good Volatility Spillover to France from Other countries")+ylim(0,1.75)

# France Spillover to other countries: bad volatility fifth

France_to_Other_bd = theta6mat[,7]
France_to_Other_bd0 = theta6mat[,9]
state = data.frame(France_to_Other_bd,France_to_Other_bd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Bad Volatility Spillover from France to Other countries")+ylim(0,1.0)

# France Spillover from other countries: bad volatility sixth

France_from_Other_bd = theta6mat[,8]
France_from_Other_bd0 = theta6mat[,10]
state = data.frame(France_from_Other_bd,France_from_Other_bd0)

df <- data.frame(x = retmat$`Date`[52:1199],
                 state)

df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color= variable)) +
  geom_line() + ggtitle("Dynamic Bad Volatility Spillover to France from Other countries")+ylim(0,1.75)


