package solver

import "errors"

// ErrorKind 区分失败原因：非法配置 / 数值失稳。
type ErrorKind int

const (
	ErrKindConfig    ErrorKind = iota // 非法配置（参数越界、未知边界类型等）
	ErrKindStability                  // 数值失稳（违反稳定性条件或产生 NaN/Inf）
)

// Error 是带类别标记的求解器错误。
type Error struct {
	Kind ErrorKind
	Msg  string
}

func (e *Error) Error() string { return e.Msg }

// KindOf 返回错误的类别；非求解器错误返回 -1。
func KindOf(err error) (ErrorKind, bool) {
	var se *Error
	if errors.As(err, &se) {
		return se.Kind, true
	}
	return 0, false
}

// IsConfigError 判定是否为非法配置错误。
func IsConfigError(err error) bool {
	k, ok := KindOf(err)
	return ok && k == ErrKindConfig
}

// IsStabilityError 判定是否为数值失稳错误。
func IsStabilityError(err error) bool {
	k, ok := KindOf(err)
	return ok && k == ErrKindStability
}

func configErr(msg string) *Error    { return &Error{Kind: ErrKindConfig, Msg: msg} }
func stabilityErr(msg string) *Error { return &Error{Kind: ErrKindStability, Msg: msg} }
